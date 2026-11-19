// Independent sparse-ledger replay. No producer code imported.
use std::collections::{HashMap,BTreeMap};
use std::io::{BufRead,BufReader,Read};
use std::fs::File;
use std::time::Instant;
type Poly=Vec<(u64,u32)>;
fn numbers(s:&str)->Vec<u64>{s.split_whitespace().map(|v|v.parse().expect("integer")).collect()}
fn deg(x:u64,r:u32,a:u32)->u32{
 let chain=x%(1u64<<r);let mut w=chain.count_ones();let mut aux=x>>r;
 for _ in 0..a{w+=(aux%64) as u32;aux/=64;}assert_eq!(aux,0,"out-of-range monomial");w
}
fn shift(src:&Poly,v:u32,m:u32,k:u32)->Poly{
 let r=m*k;let n=1u64<<k;
 src.iter().map(|&(code,c)|{
  let aux=code>>r;let mut block=code% (1u64<<r);let mut xs=Vec::new();
  for _ in 0..m{xs.push(block%n);block/=n;}
  let out=if v<r {
   let b=(v/k) as usize;xs[b]=(xs[b]+(1u64<<(v%k)))%n;
   xs.iter().enumerate().map(|(i,e)|e<<(i as u32*k)).sum::<u64>()+(aux<<r)
  }else{code+(1u64<<(r+6*(v-r)))};
  (out,c)
 }).collect()
}
fn add(h:&mut HashMap<u64,u64>,code:u64,c:u64,p:u64){
 let val=(h.get(&code).copied().unwrap_or(0)+c)%p;
 if val==0{h.remove(&code);}else{h.insert(code,val);}
}
fn closed_residual(raw:Poly,basis:&Vec<Poly>,piv:&HashMap<u64,usize>,p:u64,r:u32,a:u32)->bool{
 let mut row=BTreeMap::new();for (x,c) in raw{row.insert((deg(x,r,a),x),c as u64);}
 while let Some((&key,&c))=row.last_key_value(){
  let id=match piv.get(&key.1){Some(&id)=>id,None=>return false};
  for &(x,b) in &basis[id]{
   let key=(deg(x,r,a),x);let old=row.get(&key).copied().unwrap_or(0);
   let val=(old+p-(c*b as u64)%p)%p;
   if val==0{row.remove(&key);}else{row.insert(key,val);}
  }
 }
 true
}
fn main(){
 let args:Vec<String>=std::env::args().collect();assert!(args.len()==4||args.len()==5,"input ledger basis [closed]");
 let start=Instant::now();let mut text=String::new();File::open(&args[1]).unwrap().read_to_string(&mut text).unwrap();let data=numbers(&text);
 let p=data[0];let m=data[1] as u32;let k=data[2] as u32;let a=data[3] as u32;let d=data[4] as u32;let r=m*k;
 assert!(p>3&&p<(1<<24)&&(p-1)%(1u64<<k)==0&&r+6*a<64);
 let mut trial=2;while trial*trial<=p{assert_ne!(p%trial,0,"not prime");trial+=1;}
 let ng=data[8] as usize;let mut at=9;let mut gens=Vec::new();
 for _ in 0..ng{let n=data[at] as usize;at+=1;let mut v=Vec::new();for _ in 0..n{let x=data[at];let c=data[at+1];at+=2;assert!(c>0&&c<p);v.push((x,c as u32));}gens.push(v);}assert_eq!(at,data.len());
 let logs=BufReader::new(File::open(&args[2]).unwrap());let mut rows=BufReader::new(File::open(&args[3]).unwrap()).lines();
 let mut basis:Vec<Poly>=Vec::new();let mut piv=HashMap::new();let mut ops=0u64;
 for (id,line) in logs.lines().enumerate(){
  let line=line.unwrap();let toks:Vec<&str>=line.split_whitespace().collect();assert!(toks.len()>=6);
  assert_eq!(toks[0].parse::<usize>().unwrap(),id);let op=toks[1];let source=toks[2].parse::<usize>().unwrap();let var=toks[3].parse::<u32>().unwrap();let inv=toks[4].parse::<u64>().unwrap();let n=toks[5].parse::<usize>().unwrap();
  assert!(inv>0&&inv<p);assert_eq!(toks.len(),6+2*n);
  let raw=match op {
   "input"=>{assert!(source<ng);gens[source].clone()},
   "shift"=>{assert!(source<id&&var<r+a);assert!(basis[source].iter().all(|&(x,_)|deg(x,r,a)<d),"overdegree source");shift(&basis[source],var,m,k)},
   _=>panic!("unknown operation"),
  };
  let mut result=HashMap::new();for (x,c) in raw{assert!(deg(x,r,a)<=d);add(&mut result,x,c as u64,p);}
  for j in 0..n{let old=toks[6+2*j].parse::<usize>().unwrap();let c=toks[7+2*j].parse::<u64>().unwrap();assert!(old<id&&c>0&&c<p);
   for &(x,b) in &basis[old]{add(&mut result,x,(p-(c*b as u64)%p)%p,p);ops+=1;}
  }
  for v in result.values_mut(){*v=*v*inv%p;}
  let target=numbers(&rows.next().expect("missing row").unwrap());assert_eq!(target[0] as usize,id);let count=target[1] as usize;assert_eq!(target.len(),2+count*2);
  let mut poly=Vec::new();let mut seen=std::collections::HashSet::new();
  for j in 0..count{let x=target[2+j*2];let c=target[3+j*2];assert!(c>0&&c<p&&seen.insert(x));assert_eq!(result.get(&x).copied(),Some(c),"row identity");poly.push((x,c as u32));}
  assert_eq!(result.len(),poly.len());assert!(!poly.is_empty());
  poly.sort_by_key(|&(x,_)|std::cmp::Reverse((deg(x,r,a),x)));
  assert_eq!(poly[0].1,1);assert!(piv.insert(poly[0].0,id).is_none(),"duplicate pivot");basis.push(poly);
  if id%10000==9999{eprintln!("checked={} replay_muladds={} seconds={:.2}",id+1,ops,start.elapsed().as_secs_f64());}
 }
 assert!(rows.next().is_none(),"extra basis row");
 let check=args.len()==5;let one=piv.contains_key(&0);
 if check{
  for g in gens{if g.iter().all(|&(x,_)|deg(x,r,a)<=d){assert!(closed_residual(g,&basis,&piv,p,r,a),"missing input");}}
  for (id,row) in basis.iter().enumerate(){
   if deg(row[0].0,r,a)>=d{continue;}
   for var in 0..r+a{assert!(closed_residual(shift(row,var,m,k),&basis,&piv,p,r,a),"not closed");}
   if id%1000==0{eprintln!("closed rows checked through {} seconds={:.2}",id,start.elapsed().as_secs_f64());}
  }
 }
 println!("{{\"status\":\"passed\",\"derived_rank\":{},\"contains_one\":{},\"closed_check\":{},\"replay_muladds\":{},\"seconds\":{:.6}}}",basis.len(),one,check,ops,start.elapsed().as_secs_f64());
}
