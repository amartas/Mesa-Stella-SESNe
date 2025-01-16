#!/bin/bash

# BEGIN BLOCK
export MESA_DIR="/home/entropian/Documents/Mesa-Stella-IIP/mesa-24.08.1"
export OMP_NUM_THREADS=25
export MESASDK_ROOT="/home/entropian/mesasdk"
source "$MESASDK_ROOT/bin/mesasdk_init.sh"
export PATH="$PATH:$MESA_DIR/scripts/shmesa"
# END BLOCK

# Run make first, then run MESA
scripts=(
    mk
    rn
)

# Run
for script in "${scripts[@]}"; do
    if [ -x "$script" ]; then
        echo "Running $script..."
        bash "$script"
        if [ $? -ne 0 ]; then
            echo "Error running $script"
            exit 1
        fi
    else
        echo "Script not found or not executable: $script" # This was largely for troubleshooting
    fi
done

echo "All scripts completed"
