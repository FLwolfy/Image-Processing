import os
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
            except Exception as e:
                print(f"Access Denied: Failed to remove existing file: {target_file}, {e}")
                print(f"The file '{target_file}' is in use.")
                os.remove(module_file)
                success = False 
                
if success:
    os.rename(module_file, target_file)
    print(f"Moved {module_file} to {target_dir}")
    print(f"\nSuccess!")
else:
    print(f"\nFailed!")
