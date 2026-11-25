// Complex numbers are represented by [real, imaginary]. All angles are radians.
function unitary3(p, special=false) {
  const z=(r=0,i=0)=>[r,i], mul=([a,b],[c,d])=>[a*c-b*d,a*d+b*c], add=([a,b],[c,d])=>[a+c,b+d], conj=([a,b])=>[a,-b], phase=x=>z(Math.cos(x),Math.sin(x));
  const eye=()=>Array.from({length:3},(_,i)=>Array.from({length:3},(_,j)=>z(+(i===j))));
  const mm=(A,B)=>A.map(row=>B[0].map((_,j)=>row.reduce((s,x,k)=>add(s,mul(x,B[k][j])),z())));
  const rot=(i,j,t,f)=>{const M=eye(),c=Math.cos(t),s=Math.sin(t);M[i][i]=M[j][j]=z(c);M[i][j]=mul(z(s),phase(-f));M[j][i]=mul(z(-s),phase(f));return M;};
  const {t12=0,t13=0,t23=0,f12=0,f13=0,f23=0,a1=0,a2=0,a3=0}=p;
  const al=[a1,a2,special?-a1-a2:a3];
  const D=eye();al.forEach((x,i)=>D[i][i]=phase(x));
  const U=mm(mm(mm(D,rot(1,2,t23,f23)),rot(0,2,t13,f13)),rot(0,1,t12,f12));
  const adj=U[0].map((_,j)=>U.map(row=>conj(row[j])));
  const gram=mm(adj,U);
  const err=Math.max(...gram.flatMap((row,i)=>row.map(([r,s],j)=>Math.hypot(r-+(i===j),s))));
  const det=[[0,1,2,1],[1,2,0,1],[2,0,1,1],[2,1,0,-1],[1,0,2,-1],[0,2,1,-1]].reduce((v,[a,b,c,s])=>add(v,mul(z(s),mul(mul(U[0][a],U[1][b]),U[2][c]))),z());
  return {U,gram,err,det,al};
}
if(typeof module !== 'undefined') module.exports = { unitary3 };
