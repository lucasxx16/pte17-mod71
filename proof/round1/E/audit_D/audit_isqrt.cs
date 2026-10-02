using System;
using System.IO;
using System.Diagnostics;
public static class EAuditIsqrt {
  public static void Run(string output) {
    const long bound=18168970;
    long expected=0,maxInitial=0,maxProduct=0,corrections=0;
    var sw=Stopwatch.StartNew();
    for(long a=0;a<=bound;++a) {
      while((expected+1)*(expected+1)<=a) ++expected;
      long r=(long)Math.Sqrt((double)a);
      if(r<0 || r>4263) throw new Exception("Initial sqrt out of certified range");
      maxInitial=Math.Max(maxInitial,r);
      while(true) {
        long square=(r+1)*(r+1); maxProduct=Math.Max(maxProduct,square);
        if(square>a) break;
        ++r; ++corrections;
      }
      while(true) {
        long square=r*r; maxProduct=Math.Max(maxProduct,square);
        if(square<=a) break;
        --r; ++corrections;
      }
      if(r!=expected) throw new Exception("Incorrect integer square root");
    }
    string json="{\"status\":\"PASS\",\"all_inputs_checked\":"+(bound+1)+",\"max_input\":"+bound+",\"max_initial_r\":"+maxInitial+",\"max_correction_product\":"+maxProduct+",\"correction_steps\":"+corrections+",\"elapsed_seconds\":"+sw.Elapsed.TotalSeconds.ToString(System.Globalization.CultureInfo.InvariantCulture)+"}";
    File.WriteAllText(output,json); Console.WriteLine(json);
  }
}
