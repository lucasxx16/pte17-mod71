using System;
using System.IO;
using System.Diagnostics;
using System.Globalization;

public static class DIntegerLatticeEnumerator {
    static int N,S,R;
    static long[,] B,M;
    static long[] D,Bounds,Z,RowSums,V;
    static long[] Nodes;
    static long Leaves,OddLeaves,ExactBallLeaves,Solutions;
    static Stopwatch Clock;
    static double NextReport;
    static string DirectoryPath;
    static long FloorDiv(long a,long b) { long q=a/b, r=a%b; return q-((r<0)?1:0); }
    static long CeilDiv(long a,long b) { return -FloorDiv(-a,b); }
    static long Isqrt(long a) {
        long r=(long)Math.Sqrt((double)a);
        while ((r+1)*(r+1)<=a) ++r;
        while (r*r>a) --r;
        return r;
    }
    static void Visit(int k,long remaining,bool tailNonzero,long tailAbs,long tailParity) {
        ++Nodes[k];
        if (k==0) {
            VisitLast(remaining,tailNonzero,tailAbs,tailParity,RowSums[0]);
            return;
        }
        long center=0; for(int j=k+1;j<N;++j) center+=M[j,k]*Z[j];
        long reach=Isqrt(remaining/D[k])+tailAbs;
        long lo=Math.Max(-Bounds[k],CeilDiv(-reach-center,S));
        long hi=Math.Min(Bounds[k],FloorDiv(reach-center,S));
        if(!tailNonzero) lo=Math.Max(0,lo);
        for(long z=lo;z<=hi;++z) {
            Z[k]=z;
            long a=Math.Max(0,Math.Abs(S*z+center)-tailAbs);
            Visit(k-1,remaining-D[k]*a*a,tailNonzero||z!=0,tailAbs+Math.Abs(z),(tailParity+z*RowSums[k])&1);
        }
        if(k==N-5 && Clock.Elapsed.TotalSeconds>=NextReport) {
            Console.WriteLine("elapsed_seconds="+Clock.Elapsed.TotalSeconds.ToString("F1",CultureInfo.InvariantCulture)+" leaves="+Leaves+" nodes="+String.Join(",",Nodes));
            NextReport=Clock.Elapsed.TotalSeconds+30;
        }
    }
    static void VisitLast(long remaining,bool tailNonzero,long tailAbs,long parity,long firstsum) {
        long center=0; for(int j=1;j<N;++j) center+=M[j,0]*Z[j];
        long reach=Isqrt(remaining/D[0])+tailAbs;
        long lo=Math.Max(-Bounds[0],CeilDiv(-reach-center,S));
        long hi=Math.Min(Bounds[0],FloorDiv(reach-center,S));
        if(!tailNonzero) lo=Math.Max(0,lo);
        for(long z=lo;z<=hi;++z) {
            Z[0]=z;
            if(!tailNonzero && z==0) continue;
            ++Leaves;
            if(((parity+z*firstsum)&1)==0) continue;
            ++OddLeaves;
            long mass=0,square=0;
            bool reject=false;
            for(int i=0;i<N;++i) {
                long a=0;
                for(int j=0;j<N;++j) a+=Z[j]*B[j,i];
                V[i]=a;
                mass+=Math.Abs(a); square+=a*a;
                if(Math.Abs(a)>3 || mass>17 || square>R) { reject=true; break; }
            }
            if(reject || mass!=17) continue;
            ++ExactBallLeaves;
            ++Solutions;
            File.AppendAllText(Path.Combine(DirectoryPath,"lattice_candidate_solutions.txt"),String.Join(",",V)+Environment.NewLine);
            Console.WriteLine("SOLUTION "+String.Join(",",V));
        }
    }
    public static void Run(string directoryPath) {
        DirectoryPath=directoryPath;
        string[] tokens=File.ReadAllText(Path.Combine(directoryPath,"lattice_certificate_input.txt")).Split((char[])null,StringSplitOptions.RemoveEmptyEntries);
        int at=0;
        N=int.Parse(tokens[at++]); S=int.Parse(tokens[at++]); R=int.Parse(tokens[at++]);
        B=new long[N,N]; M=new long[N,N]; D=new long[N]; Bounds=new long[N]; Z=new long[N]; Nodes=new long[N]; RowSums=new long[N]; V=new long[N];
        for(int i=0;i<N;++i) for(int j=0;j<N;++j) B[i,j]=long.Parse(tokens[at++]);
        for(int i=0;i<N;++i) for(int j=0;j<N;++j) M[i,j]=long.Parse(tokens[at++]);
        for(int i=0;i<N;++i) D[i]=long.Parse(tokens[at++]);
        for(int i=0;i<N;++i) Bounds[i]=long.Parse(tokens[at++]);
        if(at!=tokens.Length) throw new Exception("Unexpected input tokens");
        for(int i=0;i<N;++i) for(int j=0;j<N;++j) RowSums[i]+=B[i,j];
        File.WriteAllText(Path.Combine(directoryPath,"lattice_candidate_solutions.txt"),"");
        Clock=Stopwatch.StartNew(); NextReport=30;
        Visit(N-1,(long)R*S*S*S,false,0,0);
        string result="radius_squared="+R+Environment.NewLine+"S="+S+Environment.NewLine+"integer_scaled_gram_schmidt=true"+Environment.NewLine+"one_per_sign_pair=true"+Environment.NewLine+"leaves="+Leaves+Environment.NewLine+"odd_leaves="+OddLeaves+Environment.NewLine+"solutions="+Solutions+Environment.NewLine+"nodes_by_level="+String.Join(",",Nodes)+Environment.NewLine+"elapsed_seconds="+Clock.Elapsed.TotalSeconds.ToString("F6",CultureInfo.InvariantCulture)+Environment.NewLine;
        File.WriteAllText(Path.Combine(directoryPath,"lattice_enumeration_result.txt"),result);
        Console.WriteLine(result);
    }
}
