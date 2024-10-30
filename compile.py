import os
import shutil
import glob
import subprocess

# Set up and compile the C++ extension
success = False
try:
    subprocess.run(["python", "setup.py", "build_ext", "--inplace"], check=True)
    success = True
except subprocess.CalledProcessError as e:
    print("Compile Error: ", e)

# Move compiled module files to the modules directory
if success:
    target_dir = 'modules'
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)

    module_files = glob.glob("*.so") + glob.glob("*.pyd")
    for module_file in module_files:
        target_file = os.path.join(target_dir, os.path.basename(module_file))
        
        if os.path.exists(target_file):
            try:
                os.remove(target_file)
                print(f"Removed existing file: {target_file}")
                os.rename(module_file, target_file)
                print(f"Successfully compiled and moved {target_file}")
            except Exception as e:
                print(f"Failed to remove existing file: {target_file}\nError: {e}")
                os.remove(module_file)
