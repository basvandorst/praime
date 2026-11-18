const fs=require("fs");
const html=fs.readFileSync(process.argv[2],"utf8");
const grab=name=>{const i=html.indexOf("function "+name+"(");if(i<0)throw new Error(name);const j=html.indexOf("\n",i);return html.slice(i,j);};
let src="";for(const f of["chJac","chIsqrt","duoPow"])src+=grab(f)+"\n";
const i=html.indexOf("function n2Start("),j=html.indexOf("const ECC=new Map();",i);src+=html.slice(i,j);
const api=new Function(src+";return{n2Test,n2Start};")();
const out=[];for(const line of fs.readFileSync(process.argv[3],"utf8").split("\n")){const s=line.trim();if(!s)continue;const n=BigInt(s),Q=api.n2Test(n);out.push([s,Q.a,Q.c,Q.cells.length?Q.cells[0].toString():"",Q.cells.length?Q.cells[1].toString():"",Q.ok?1:0,Q.msg]);}
fs.writeFileSync(process.argv[4],JSON.stringify(out));
