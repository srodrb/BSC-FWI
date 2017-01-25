version 0
	- compilador intel 15 o superior
	- compilador gcc 6.2.0

para la version del offload
	- module load ompss/stable
	- module load intel/16.0.0
	- module load gcc/6.2.0

para ejecutar en KNLs
        - fichero slurm.sh reserva una hora de nodo
        - una vez que asigna un nodo (squeue)
        - hacer ssh al nodo asignado (ssh knl0x, x puede ser 1 ,2 ,3, 4)
        - ya disponible para lanzar
