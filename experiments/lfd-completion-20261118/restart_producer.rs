// Sparse immutable echelon closure with a replayable derivation for every pivot.
// Input format: p m k auxiliaries D max_seconds max_rows max_terms
// followed by generator count and (term count, code coefficient pairs) per input.
use std::collections::{HashMap,VecDeque};
use std::io::{self,Read,Write,BufWriter,BufReader,BufRead};
use std::fs::File;
use std::time::Instant;
type Row=Vec<(u32,u32)>;
// A gzip child receives only proof bytes; close and wait before reporting success.
struct ProofOutput { writer:Option<BufWriter<Box<dyn Write>>>, child:Option<std::process::Child> }
impl ProofOutput {
 fn create(path:&str)->Self {
  let file=File::create(path).unwrap();
  if path.ends_with(".gz") {
   let mut child=std::process::Command::new("gzip").args(["-n","-1","-c"])
    .stdin(std::process::Stdio::piped()).stdout(file).spawn().unwrap();
   let input=child.stdin.take().unwrap();
   Self{writer:Some(BufWriter::new(Box::new(input))),child:Some(child)}
  } else {Self{writer:Some(BufWriter::new(Box::new(file))),child:None}}
 }
 fn finish(&mut self) {
  self.flush().unwrap();drop(self.writer.take());
  if let Some(mut child)=self.child.take(){assert!(child.wait().unwrap().success(),"gzip failed");}
 }
}
impl Write for ProofOutput {
 fn write(&mut self,b:&[u8])->io::Result<usize>{self.writer.as_mut().unwrap().write(b)}
 fn flush(&mut self)->io::Result<()>{self.writer.as_mut().unwrap().flush()}
}

// Read compressed proof prefixes without expanding them on disk.
fn read_proof(path:&str)->(Box<dyn Read>,Option<std::process::Child>){
 if path.ends_with(".gz"){
  let mut child=std::process::Command::new("gzip").args(["-d","-c",path]).stdout(std::process::Stdio::piped()).spawn().unwrap();
  let pipe=child.stdout.take().unwrap();(Box::new(pipe),Some(child))
 }else{(Box::new(File::open(path).unwrap()),None)}
}
fn finish_read(child:Option<std::process::Child>){if let Some(mut c)=child{assert!(c.wait().unwrap().success(),"seed gzip decode failed");}}
fn weight(x:u64,r:u32,a:u32)->u32 {
 let mut w=(x&((1u64<<r)-1)).count_ones();
 for j in 0..a {w+=((x>>(r+6*j))&63) as u32;} w
}
fn power(mut a:u64,mut n:u64,p:u64)->u32 {let mut z=1;while n>0 {if n&1!=0 {z=z*a%p;}a=a*a%p;n>>=1;} z as u32}
fn shifted(x:u64,v:u32,m:u32,k:u32)->u64 {
 let r=m*k;
 if v<r {let off=(v/k)*k;let mask=((1u64<<k)-1)<<off;
 (x&!mask)|(((((x>>off)&((1u64<<k)-1))+(1u64<<(v%k)))&((1u64<<k)-1))<<off)
 }else {x+(1u64<<(r+6*(v-r)))}
}
const NONE:usize=usize::MAX;
fn choose(n:u32,k:u32)->u128 {if k>n{return 0;}let mut z=1;for j in 0..k{z=z*(n-j) as u128/(j+1) as u128;}z}
fn aux(codes:&mut Vec<u64>,x:u64,j:u32,a:u32,r:u32,left:u32){
 if j==a{codes.push(x);return;}
 for e in 0..=left {aux(codes,x|((e as u64)<<(r+6*j)),j+1,a,r,left-e);}
}
fn chains(codes:&mut Vec<u64>,x:u64,next:u32,r:u32,a:u32,left:u32){
 aux(codes,x,0,a,r,left);
 if left>0 {for bit in next..r{chains(codes,x|(1u64<<bit),bit+1,r,a,left-1);}}
}
struct Space {codes:Vec<u64>,index:HashMap<u64,u32>,weights:Vec<u32>,piv:Vec<usize>,values:Vec<u32>,bits:Vec<u64>}
impl Space {
 fn new(r:u32,a:u32,d:u32)->Self {
  let count:u128=(0..=d.min(r)).map(|w|choose(r,w)*choose(d-w+a,a)).sum();
  assert!(count<=2_000_000,"indexed-space cap exceeded: {} columns",count);
  let mut codes=Vec::with_capacity(count as usize);chains(&mut codes,0,0,r,a,d);
  codes.sort_unstable_by_key(|&x|std::cmp::Reverse((weight(x,r,a),x)));
  assert_eq!(codes.len(),count as usize);let n=codes.len();
  let index=codes.iter().enumerate().map(|(i,&x)|(x,i as u32)).collect();
  let weights=codes.iter().map(|&x|weight(x,r,a)).collect();
  Self{codes,index,weights,piv:vec![NONE;n],values:vec![0;n],bits:vec![0;(n+63)/64]}
 }
 fn add(&mut self,x:u32,c:u32,p:u32){
  let j=x as usize;let z=((self.values[j] as u64+c as u64)%p as u64) as u32;
  self.values[j]=z;let bit=1u64<<(j%64);
  if z==0{self.bits[j/64]&=!bit;}else{self.bits[j/64]|=bit;}
 }
 fn leading(&self,cursor:&mut usize)->Option<usize>{
  while *cursor<self.bits.len(){let b=self.bits[*cursor];if b!=0{return Some(*cursor*64+b.trailing_zeros() as usize);}*cursor+=1;}None
 }
}
fn absorb(raw:Row,source:(&str,usize,u32),p:u32,d:u32,space:&mut Space,
 basis:&mut Vec<Row>,queue:&mut VecDeque<usize>,log:&mut ProofOutput,out:&mut ProofOutput)->(bool,u64){
 for (x,c) in raw{space.add(x,c,p);}
 let mut reductions=Vec::new();let mut ops=0;let mut cursor=0;
 while let Some(x)=space.leading(&mut cursor){
  let i=space.piv[x];if i==NONE{break;}
  let c=space.values[x];reductions.push((i,c));
  for &(y,b) in &basis[i]{let z=(c as u64*b as u64%p as u64) as u32;space.add(y,if z==0{0}else{p-z},p);ops+=1;}
 }
 let x=match space.leading(&mut cursor){Some(x)=>x,None=>return (false,ops)};
 let inv=power(space.values[x] as u64,p as u64-2,p as u64);let mut new=Vec::new();
 for word in cursor..space.bits.len(){let mut b=space.bits[word];
  while b!=0{let j=word*64+b.trailing_zeros() as usize;let v=space.values[j];
   new.push((j as u32,(v as u64*inv as u64%p as u64) as u32));space.values[j]=0;b&=b-1;}
  space.bits[word]=0;
 }
 let id=basis.len();space.piv[x]=id;if space.weights[x]<d{queue.push_back(id);}
 write!(log,"{} {} {} {} {} {}",id,source.0,source.1,source.2,inv,reductions.len()).unwrap();
 for (i,c) in reductions{write!(log," {} {}",i,c).unwrap();}writeln!(log).unwrap();
 write!(out,"{} {}",id,new.len()).unwrap();for &(j,c) in &new{write!(out," {} {}",space.codes[j as usize],c).unwrap();}writeln!(out).unwrap();
 basis.push(new);(space.codes[x]==0,ops)
}
fn main(){
 let args:Vec<String>=std::env::args().collect();assert!(args.len()==4||args.len()==6||args.len()==7,"input ledger basis [seed_ledger seed_basis [schedule_products]]");
 let mut text=String::new();File::open(&args[1]).unwrap().read_to_string(&mut text).unwrap();let mut it=text.split_whitespace();
 macro_rules! n {()=>{it.next().unwrap().parse::<u64>().unwrap()};}
 let p=n!() as u32;let m=n!() as u32;let k=n!() as u32;let a=n!() as u32;let d=n!() as u32;
 let seconds=n!();let max_rows=n!() as usize;let max_terms=n!();let r=m*k;
 assert!(p>3&&p<1<<24&&r>0&&r+6*a<64&&d<63&&k>0);
 let ng=n!();let mut gens=Vec::new();for _ in 0..ng{let nt=n!();let mut row=Vec::new();for _ in 0..nt{let x=n!();let c=n!();assert!(c<p as u64);row.push((x,c as u32));}gens.push(row);}assert!(it.next().is_none());
 // Preflight before creating proof files; no implicit truncated monomial domain.
 let mut space=Space::new(r,a,d);let start=Instant::now();
 let mut log=ProofOutput::create(&args[2]);let mut out=ProofOutput::create(&args[3]);
 let mut basis=Vec::new();let mut queue=VecDeque::new();let mut one=false;let mut ops=0;let mut stored=0u64;
 if args.len()>=6 {
  // Seed membership must be checked by the caller; the combined proof is
  // independently replayed again. Requeue every admissible seed row so
  // restart cannot silently skip products absent from the retained log.
  let (reader,child)=read_proof(&args[5]);
  for line in BufReader::new(reader).lines(){let line=line.unwrap();
   let v:Vec<u64>=line.split_whitespace().map(|x|x.parse().unwrap()).collect();
   assert!(v.len()>=4&&v[0] as usize==basis.len()&&v.len()==2+2*v[1] as usize);
   let mut row:Row=Vec::new();
   for pair in v[2..].chunks_exact(2){let j=*space.index.get(&pair[0]).expect("seed outside domain");
    assert!(pair[1]>0&&pair[1]<p as u64);
    if let Some(&(old,_))=row.last(){assert!(old<j,"seed ordering");}
    row.push((j,pair[1] as u32));
   }
   assert_eq!(row[0].1,1);let pivot=row[0].0 as usize;
   assert_eq!(space.piv[pivot],NONE,"duplicate seed pivot");space.piv[pivot]=basis.len();
   if space.weights[pivot]<d{queue.push_back(basis.len());}
   one|=space.codes[pivot]==0;stored+=row.len() as u64;
   assert!(basis.len()<max_rows&&stored<=max_terms,"seed exceeds resource budget");
   basis.push(row);writeln!(out,"{}",line).unwrap();
  }finish_read(child);
  let (mut reader,child)=read_proof(&args[4]);io::copy(&mut reader,&mut log).unwrap();finish_read(child);
 }
 let seed_rank=basis.len();
 for (i,g) in gens.into_iter().enumerate(){if g.iter().any(|&(x,_)|weight(x,r,a)>d){continue;}
  let raw=g.into_iter().map(|(x,c)|(*space.index.get(&x).expect("input monomial outside domain"),c)).collect();
  let (hit,o)=absorb(raw,("input",i,0),p,d,&mut space,&mut basis,&mut queue,&mut log,&mut out);one|=hit;ops+=o;
 }
 let mut first_var=0;
 if args.len()==7 {
  let done=args[6].parse::<usize>().unwrap();let nvars=(r+a) as usize;
  let complete=done/nvars;first_var=(done%nvars) as u32;
  assert!(complete<=queue.len(),"seed cursor beyond queue");
  for _ in 0..complete{queue.pop_front().unwrap();}
  assert!(first_var==0||!queue.is_empty(),"partial cursor beyond queue");
 }
 let mut products=0;let mut status="closed";
 'outer:while !one{let id=match queue.pop_front(){Some(i)=>i,None=>break};
  let begin=first_var;first_var=0;
  for v in begin..r+a{
   if start.elapsed().as_secs()>=seconds||basis.len()>=max_rows||stored>=max_terms{status="capped";break 'outer;}
   let old=basis.len();let raw=basis[id].iter().map(|&(j,c)|{
    let code=shifted(space.codes[j as usize],v,m,k);(*space.index.get(&code).expect("shift outside domain"),c)
   }).collect();
   let (hit,o)=absorb(raw,("shift",id,v),p,d,&mut space,&mut basis,&mut queue,&mut log,&mut out);one|=hit;ops+=o;products+=1;
   if basis.len()>old{stored+=basis.last().unwrap().len() as u64;}if one{break;}
   if products%10000==0{eprintln!("products={} rank={} terms={} queued={} seconds={:.1}",products,basis.len(),stored,queue.len(),start.elapsed().as_secs_f64());log.flush().unwrap();out.flush().unwrap();}
  }
 }
 if one{status="unit-witness";}log.finish();out.finish();
 eprintln!("seed_rank={} final_rank={}",seed_rank,basis.len());
 println!("{{\"status\":\"{}\",\"p\":{},\"m\":{},\"k\":{},\"D\":{},\"derived_rank\":{},\"contains_one\":{},\"products\":{},\"terms\":{},\"muladds\":{},\"seconds\":{:.6}}}",status,p,m,k,d,basis.len(),one,products,stored,ops,start.elapsed().as_secs_f64());
 io::stdout().flush().unwrap();
}
