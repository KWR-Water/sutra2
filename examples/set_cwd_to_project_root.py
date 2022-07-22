import os
from pathlib import Path
exists = False
n=0
while (not exists) :
    try:
        import research
        exists = True
    except ModuleNotFoundError:
        # traverse one level up, and see if it can be imported from there.
        os.chdir('..')

    if n==8:
        raise ModuleNotFoundError(f"Could not import from project root, tried {n} directories levels above current file's directory")

    n += 1

project_root = Path(os.getcwd())