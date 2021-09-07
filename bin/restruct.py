from __future__ import division
import shlex, subprocess
import commands

cmd="rm -rf res/*"
sta,out=commands.getstatusoutput(cmd)
cmd="rm -rf HeteroMirPred/res_frag_*"
sta,out=commands.getstatusoutput(cmd)

f=open("record_log.txt","w")
f.writelines(["100001\n"])
f.close()