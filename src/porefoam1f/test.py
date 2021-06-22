
import os; ''' ========== set up paths  =========== '''
if not ("msRoot" in os.environ): 
  print("try again after running:\nsource .../src/bashrc"); exit(-1);
from msrc import *;   
DbgMsg('============ Ignore above messages ===============')









with open("voxcylinder.mhd", 'w') as f1:
	f1.write("""DimSize = 20 20 20 
				Offset =      0    0    0
				replaceRange 0 255 1
				reset  dx 1 1 1
				Paint cylinder 0 10 10   20 10 10  5
				reset  dx 1e-6 1e-6 1e-6
				""");#ElementDataFile = NO_READ

runSh('.', "rm -r voxcylinder/*");
runSh('.', "AllRunImagePar voxcylinder.mhd");
exit(fileFloatDiffersFrom("voxcylinder/voxcylinder-1-X/summary_voxcylinder-1-X.txt","K_x= ",2.44436e-12))

          
