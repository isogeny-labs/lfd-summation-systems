// Bounded divisor reduction with a membership derivation for every retained row.
// Empty critical-pair queue is provisional until the bounded criterion is proved.
use std::collections::{BTreeMap,BinaryHeap,HashMap};
use std::cmp::Reverse;
use std::fs::File;
use std::io::{Read,Write,BufWriter};
use std::time::Instant;
type Row=Vec<(u64,u32)>;
type Acc=BTreeMap<(u32,u64),u32>;
#[derive(Clone,Copy)]
struct Ring{p:u32,m:u32,k:u32,a:u32,d:u32,r:u32}
impl Ring{
 fn degree(&self,x:u64)->u32{let mut w=(x&((1u64<<self.r)-1)).count_ones();for j in 0..self.a{w+=((x>>(self.r+6*j))&63) as u32;}w}
 fn mul(&self,x:u64,y:u64)->u64{let mask=(1u64<<self.k)-1;let mut z=0;for b in 0..self.m{let shift=b*self.k;z|=((((x>>shift)&mask)+((y>>shift)&mask))&mask)<<shift;}for j in 0..self.a{let shift=self.r+6*j;let e=((x>>shift)&63)+((y>>shift)&63);assert!(e<64);z|=e<<shift;}z}
 fn quotient(&self,x:u64,y:u64)->Option<u64>{let mask=(1u64<<self.r)-1;if (x&y&mask)!=(y&mask){return None;}let mut z=(x^y)&mask;for j in 0..self.a{let shift=self.r+6*j;let ex=(x>>shift)&63;let ey=(y>>shift)&63;if ex<ey{return None;}z|=(ex-ey)<<shift;}Some(z)}
 fn lcm(&self,x:u64,y:u64)->u64{let mut z=(x|y)&((1u64<<self.r)-1);for j in 0..self.a{let shift=self.r+6*j;z|=((x>>shift)&63).max((y>>shift)&63)<<shift;}z}
 fn coprime(&self,x:u64,y:u64)->bool{if (x&y&((1u64<<self.r)-1))!=0{return false;}for j in 0..self.a{let shift=self.r+6*j;if ((x>>shift)&63)>0&&((y>>shift)&63)>0{return false;}}true}
 fn aux_search(&self,x:u64,prefix:u64,j:u32,lookup:&HashMap<u64,usize>,best:&mut Option<usize>){
  if *best==Some(0){return;}
  if j==self.a{if let Some(&id)=lookup.get(&prefix){if best.map_or(true,|old|id<old){*best=Some(id);}}return;}
  let shift=self.r+6*j;for e in 0..=((x>>shift)&63){self.aux_search(x,prefix|(e<<shift),j+1,lookup,best);}
 }
 fn divisor(&self,x:u64,lookup:&HashMap<u64,usize>)->Option<usize>{let mask=x&((1u64<<self.r)-1);let mut sub=mask;let mut best=None;loop{self.aux_search(x,sub,0,lookup,&mut best);if sub==0||best==Some(0){break;}sub=(sub-1)&mask;}best}
 fn put(&self,row:&mut Acc,x:u64,c:u32){if c==0{return;}let key=(self.degree(x),x);let value=((row.get(&key).copied().unwrap_or(0) as u64+c as u64)%self.p as u64) as u32;if value==0{row.remove(&key);}else{row.insert(key,value);}}
 fn shifted(&self,row:&Row,h:u64)->Row{assert!(self.degree(row[0].0)+self.degree(h)<=self.d);row.iter().map(|&(x,c)|(self.mul(x,h),c)).collect()}
}
fn inverse(mut a:u64,p:u64)->u32{let mut n=p-2;let mut z=1;while n!=0{if n&1!=0{z=z*a%p;}a=a*a%p;n>>=1;}z as u32}
#[derive(Clone,Copy)]
enum Source{Input(usize),Carry(usize,u32),Pair(usize,usize)}
struct Engine{
 ring:Ring,basis:Vec<Row>,leads:Vec<u64>,lookup:HashMap<u64,usize>,queue:BinaryHeap<Reverse<(u32,u64,u8,usize,usize)>>,
 serial:u64,processed:u64,ops:u64,terms:u64,prune_coprime:bool,frontier:bool,skipped:u64,max_rows:usize,max_terms:u64,max_queue:usize,seconds:u64,start:Instant,log:BufWriter<File>,nodes:usize,
}
impl Engine{
 fn budget(&self)->Result<(),()>{if self.start.elapsed().as_secs()>=self.seconds||self.basis.len()>=self.max_rows||self.terms>=self.max_terms||self.queue.len()>=self.max_queue{Err(())}else{Ok(())}}
 fn enqueue(&mut self,degree:u32,kind:u8,i:usize,j:usize)->Result<(),()>{self.budget()?;self.queue.push(Reverse((if self.frontier{self.ring.d-degree}else{degree},if self.frontier{u64::MAX-self.serial}else{self.serial},kind,i,j)));self.serial+=1;Ok(())}
 fn emit(&mut self,source:Source,inv:u32,ops:&[(usize,u64,u32)],row:&Row){
  if self.nodes>0{write!(self.log,",").unwrap();}
  match source{Source::Input(i)=>write!(self.log,"{{\"source\":[\"input\",{}]",i),Source::Carry(i,v)=>write!(self.log,"{{\"source\":[\"carry\",{},{}]",i,v),Source::Pair(i,j)=>write!(self.log,"{{\"source\":[\"pair\",{},{}]",i,j)}.unwrap();
  write!(self.log,",\"inverse\":{},\"reductions\":[",inv).unwrap();for (n,(i,h,c)) in ops.iter().enumerate(){if n>0{write!(self.log,",").unwrap();}write!(self.log,"[{},{},{}]",i,h,c).unwrap();}
  write!(self.log,"],\"row\":[").unwrap();let mut sorted=row.clone();sorted.sort_unstable_by_key(|&(x,_)|x);for(n,(x,c))in sorted.iter().enumerate(){if n>0{write!(self.log,",").unwrap();}write!(self.log,"[{},{}]",x,c).unwrap();}write!(self.log,"]}}").unwrap();self.nodes+=1;
 }
 fn add(&mut self,raw:Row,source:Source)->Result<(),()>{
  let ring=self.ring;let mut row=Acc::new();for(x,c)in raw{assert!(ring.degree(x)<=ring.d);ring.put(&mut row,x,c);}
  let mut rem=Vec::new();let mut operations=Vec::new();
  while let Some((&key,&c))=row.last_key_value(){
   self.budget()?;let x=key.1;
   if let Some(j)=ring.divisor(x,&self.lookup){
    let h=ring.quotient(x,self.leads[j]).unwrap();assert!(ring.degree(self.leads[j])+ring.degree(h)<=ring.d);operations.push((j,h,c));
    for (index,&(y,b)) in self.basis[j].iter().enumerate(){let z=ring.mul(y,h);if index==0{assert_eq!(z,x);assert_eq!(b,1);}else{assert!((ring.degree(z),z)<key,"leading multiplication");}
     let v=(c as u64*b as u64%ring.p as u64) as u32;ring.put(&mut row,z,if v==0{0}else{ring.p-v});self.ops+=1;if self.ops%4096==0{self.budget()?;}
    }
   }else{row.remove(&key);rem.push((x,c));}
  }
  if rem.is_empty(){return Ok(());}
  let lead=rem[0].0;let inv=inverse(rem[0].1 as u64,ring.p as u64);for(_,c)in &mut rem{*c=(*c as u64*inv as u64%ring.p as u64) as u32;}
  self.emit(source,inv,&operations,&rem);let i=self.basis.len();self.terms+=rem.len() as u64;self.basis.push(rem);self.leads.push(lead);assert!(self.lookup.insert(lead,i).is_none());
  if lead==0{return Ok(());}
  let degree=ring.degree(lead);
  if degree<ring.d{for v in 0..ring.r{if lead>>v&1!=0{self.enqueue(degree+1,0,i,v as usize)?;}}}
  for j in 0..i{let deg=ring.degree(ring.lcm(lead,self.leads[j]));if deg<=ring.d{if self.prune_coprime&&ring.coprime(lead,self.leads[j]){self.skipped+=1;}else{self.enqueue(deg,1,i,j)?;}}if j%4096==0{self.budget()?;}}
  Ok(())
 }
 fn run(&mut self,gens:Vec<Row>)->Result<(),()>{
  for(i,g)in gens.into_iter().enumerate(){if g.iter().all(|&(x,_)|self.ring.degree(x)<=self.ring.d){self.add(g,Source::Input(i))?;}if self.lookup.contains_key(&0){return Ok(());}}
  while let Some(Reverse((_,_,kind,i,j)))=self.queue.pop(){
   self.budget()?;let raw;if kind==0{raw=self.ring.shifted(&self.basis[i],1u64<<j);}else{
    let lcm=self.ring.lcm(self.leads[i],self.leads[j]);raw={let mut left=self.ring.shifted(&self.basis[i],self.ring.quotient(lcm,self.leads[i]).unwrap());let right=self.ring.shifted(&self.basis[j],self.ring.quotient(lcm,self.leads[j]).unwrap());left.extend(right.into_iter().map(|(x,c)|(x,self.ring.p-c)));left};
   }
   self.add(raw,if kind==0{Source::Carry(i,j as u32)}else{Source::Pair(i,j)})?;self.processed+=1;
   if self.lookup.contains_key(&0){return Ok(());}
   if self.processed%1000==0{eprintln!("pairs={} compact_rows={} terms={} queue={} seconds={:.2}",self.processed,self.basis.len(),self.terms,self.queue.len(),self.start.elapsed().as_secs_f64());self.log.flush().unwrap();}
  }Ok(())
 }
}
fn main(){
 let args:Vec<String>=std::env::args().collect();assert!(args.len()>=3&&args.len()<=6,"input proof [max_queue [coprime] [frontier]]");for flag in args.iter().skip(4){assert!(flag=="coprime"||flag=="frontier");};let mut text=String::new();File::open(&args[1]).unwrap().read_to_string(&mut text).unwrap();let mut it=text.split_whitespace();
 macro_rules! n{()=>{it.next().unwrap().parse::<u64>().unwrap()};}
 let p=n!() as u32;let m=n!() as u32;let k=n!() as u32;let a=n!() as u32;let d=n!() as u32;let seconds=n!();let max_rows=n!() as usize;let max_terms=n!();let r=m*k;
 assert!(p>3&&p<(1<<24)&&k>0&&r>0&&r+6*a<64&&d<63);let ring=Ring{p,m,k,a,d,r};
 let count=n!();let mut gens=Vec::new();for _ in 0..count{let n=n!();let mut g=Vec::new();for _ in 0..n{let x=n!();let c=n!();assert!(x>>(r+6*a)==0&&c<p as u64);g.push((x,c as u32));}gens.push(g);}assert!(it.next().is_none());
 let mut engine=Engine{ring,basis:Vec::new(),leads:Vec::new(),lookup:HashMap::new(),queue:BinaryHeap::new(),serial:0,processed:0,ops:0,terms:0,prune_coprime:args.iter().skip(4).any(|s|s=="coprime"),frontier:args.iter().skip(4).any(|s|s=="frontier"),skipped:0,max_rows,max_terms,max_queue:if args.len()>=4{args[3].parse().unwrap()}else{2_000_000},seconds,start:Instant::now(),log:BufWriter::new(File::create(&args[2]).unwrap()),nodes:0};
 write!(engine.log,"{{\"trace\":[").unwrap();let capped=engine.run(gens).is_err();write!(engine.log,"]}}\n").unwrap();engine.log.flush().unwrap();let one=engine.lookup.contains_key(&0);
 println!("{{\"status\":\"{}\",\"contains_one\":{},\"compact_rows\":{},\"terms\":{},\"pairs_processed\":{},\"queued\":{},\"coprime_pairs_skipped\":{},\"muladds\":{},\"seconds\":{:.6},\"claim\":\"membership pending replay; compact closedness criterion unproved\"}}",if one{"unit-witness"}else if capped{"capped"}else{"provisional-complete"},one,engine.basis.len(),engine.terms,engine.processed,engine.queue.len(),engine.skipped,engine.ops,engine.start.elapsed().as_secs_f64());
}
