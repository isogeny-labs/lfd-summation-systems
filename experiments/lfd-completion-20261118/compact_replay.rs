// Independent compact proof replay. No producer algebra is imported.
use std::collections::{HashMap,HashSet};
use std::fs::File;
use std::io::{Read,BufRead,BufReader};
use std::time::Instant;
type Poly=Vec<(u32,u32)>;
fn nums(s:&str)->Vec<u64>{s.split_whitespace().map(|v|v.parse().expect("integer token")).collect()}
struct Columns{map:HashMap<u64,u32>,codes:Vec<u64>,values:Vec<u32>,marked:Vec<bool>,dirty:Vec<usize>,nonzero:usize}
impl Columns{
 fn new()->Self{Self{map:HashMap::new(),codes:Vec::new(),values:Vec::new(),marked:Vec::new(),dirty:Vec::new(),nonzero:0}}
 fn index(&mut self,x:u64)->u32{if let Some(&j)=self.map.get(&x){return j;}assert!(self.codes.len()<u32::MAX as usize);let j=self.codes.len() as u32;self.map.insert(x,j);self.codes.push(x);self.values.push(0);self.marked.push(false);j}
 fn add(&mut self,j:u32,c:u32,p:u64){let j=j as usize;let old=self.values[j];let value=((old as u64+c as u64)%p) as u32;if old==0&&value!=0{self.nonzero+=1;}if old!=0&&value==0{self.nonzero-=1;}self.values[j]=value;if !self.marked[j]{self.marked[j]=true;self.dirty.push(j);}}
 fn clear(&mut self){for j in self.dirty.drain(..){self.values[j]=0;self.marked[j]=false;}self.nonzero=0;}
}
#[derive(Clone,Copy)]
struct Algebra{p:u64,m:u32,k:u32,a:u32,d:u32,r:u32}
impl Algebra{
 fn coordinates(&self,x:u64)->(Vec<u64>,Vec<u64>){assert_eq!(x>>(self.r+6*self.a),0,"monomial range");let n=1u64<<self.k;let chain=(0..self.m).map(|b|(x>>(b*self.k))%n).collect();let aux=(0..self.a).map(|j|(x>>(self.r+6*j))%64).collect();(chain,aux)}
 fn degree(&self,x:u64)->u32{let (chain,aux)=self.coordinates(x);chain.iter().map(|v|v.count_ones()).sum::<u32>()+aux.iter().sum::<u64>() as u32}
 fn product(&self,x:u64,h:u64)->u64{let (xs,es)=self.coordinates(x);let (ys,fs)=self.coordinates(h);let mut out=0;for b in 0..self.m as usize{out+=((xs[b]+ys[b])%(1u64<<self.k))<<(b as u32*self.k);}for j in 0..self.a as usize{let e=es[j]+fs[j];assert!(e<64,"auxiliary overflow");out+=e<<(self.r+6*j as u32);}out}
 fn pair_multipliers(&self,x:u64,y:u64)->(u64,u64){
  let (xs,ex)=self.coordinates(x);let (ys,ey)=self.coordinates(y);let mut hx=0;let mut hy=0;let mut degree=0;
  for b in 0..self.m as usize{let union=xs[b]|ys[b];degree+=union.count_ones();hx+=(union-xs[b])<<(b as u32*self.k);hy+=(union-ys[b])<<(b as u32*self.k);}
  for j in 0..self.a as usize{let e=ex[j].max(ey[j]);degree+=e as u32;hx+=(e-ex[j])<<(self.r+6*j as u32);hy+=(e-ey[j])<<(self.r+6*j as u32);}
  assert!(degree<=self.d,"pair degree budget");(hx,hy)
 }
 fn add_multiple(&self,row:&Poly,lead:u64,h:u64,c:u64,columns:&mut Columns,ops:&mut u64){
  assert!(self.degree(lead)+self.degree(h)<=self.d,"product degree budget");
  for &(j,b) in row{
   let dest=if h==0{j}else{let code=self.product(columns.codes[j as usize],h);assert!(self.degree(code)<=self.d);columns.index(code)};
   columns.add(dest,(c*b as u64%self.p) as u32,self.p);*ops+=1;
  }
 }
}
fn main(){
 let args:Vec<String>=std::env::args().collect();assert_eq!(args.len(),3,"input text-proof");let start=Instant::now();let mut input=String::new();File::open(&args[1]).unwrap().read_to_string(&mut input).unwrap();let data=nums(&input);assert!(data.len()>=9);
 let p=data[0];let m=data[1] as u32;let k=data[2] as u32;let a=data[3] as u32;let d=data[4] as u32;let r=m*k;
 assert!(p>3&&p<(1<<24)&&m>0&&k>0&&r+6*a<64&&d<63&&(p-1)%(1u64<<k)==0);
 let mut trial=2;while trial*trial<=p{assert_ne!(p%trial,0,"prime field required");trial+=1;}
 let algebra=Algebra{p,m,k,a,d,r};let ng=data[8] as usize;let mut at=9;let mut gens=Vec::new();
 for _ in 0..ng{let n=data[at] as usize;at+=1;let mut row=Vec::new();let mut seen=HashSet::new();for _ in 0..n{let x=data[at];let c=data[at+1];at+=2;assert!(c>0&&c<p&&seen.insert(x));algebra.coordinates(x);row.push((x,c as u32));}gens.push(row);}assert_eq!(at,data.len());
 let mut lines=BufReader::new(File::open(&args[2]).unwrap()).lines();let count=lines.next().expect("proof count").unwrap().parse::<usize>().unwrap();let mut basis:Vec<Poly>=Vec::new();let mut leads=Vec::new();let mut columns=Columns::new();let mut ops=0;let mut one=false;
 for id in 0..count{
  let values=nums(&lines.next().expect("missing node").unwrap());assert!(values.len()>=7);assert_eq!(values[0] as usize,id);let kind=values[1];let i=values[2] as usize;let j=values[3] as usize;let inv=values[4];let reductions=values[5] as usize;assert!(inv>0&&inv<p);
  match kind{
   0=>{assert!(i<ng&&j==0);for &(x,c) in &gens[i]{assert!(algebra.degree(x)<=d,"input degree");let dest=columns.index(x);columns.add(dest,c,p);}},
   1=>{assert!(i<id&&j<r as usize,"carry reference");algebra.add_multiple(&basis[i],leads[i],1u64<<j,1,&mut columns,&mut ops);},
   2=>{assert!(i<id&&j<id,"pair reference");let (hi,hj)=algebra.pair_multipliers(leads[i],leads[j]);algebra.add_multiple(&basis[i],leads[i],hi,1,&mut columns,&mut ops);algebra.add_multiple(&basis[j],leads[j],hj,p-1,&mut columns,&mut ops);},
   _=>panic!("source kind")
  }
  let mut at=6;for _ in 0..reductions{assert!(at+3<=values.len());let old=values[at] as usize;let h=values[at+1];let c=values[at+2];at+=3;assert!(old<id&&c>0&&c<p,"reducer reference");algebra.add_multiple(&basis[old],leads[old],h,p-c,&mut columns,&mut ops);}
  assert!(at<values.len());let terms=values[at] as usize;at+=1;assert_eq!(values.len(),at+2*terms);assert!(terms>0);let mut target=Vec::new();let mut seen=HashSet::new();let mut lead=None;
  for _ in 0..terms{let x=values[at];let c=values[at+1];at+=2;let degree=algebra.degree(x);assert!(c>0&&c<p&&degree<=d&&seen.insert(x));let dest=columns.index(x);assert_eq!(columns.values[dest as usize] as u64*inv%p,c,"row identity");target.push((dest,c as u32));if lead.map_or(true,|(w,y,_)|(degree,x)>(w,y)){lead=Some((degree,x,c));}}
  assert_eq!(columns.nonzero,terms,"extra residual term");let (_,lead_code,c)=lead.unwrap();assert_eq!(c,1,"monic output");if lead_code==0{one=true;}
  basis.push(target);leads.push(lead_code);columns.clear();if id%100==99{eprintln!("checked={} muladds={} seconds={:.2}",id+1,ops,start.elapsed().as_secs_f64());}
 }
 assert!(lines.next().is_none(),"extra proof node");println!("{{\"status\":\"passed\",\"verified_rows\":{},\"contains_one\":{},\"closed_check\":false,\"replay_muladds\":{},\"seconds\":{:.6},\"scope\":\"independent compact row membership only\"}}",basis.len(),one,ops,start.elapsed().as_secs_f64());
}
