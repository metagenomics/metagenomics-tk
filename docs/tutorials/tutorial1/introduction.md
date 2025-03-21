



Todos:
* Introduce the Toolkti

* Introduce the data to process
* Download gtdb!
* Set a working directory to a volume : ~/mgcourse/ should link to /vol/mgcourse !!!

* Explain sections of the config

* Explain full pipeline and module specific executions.

* Explain the usual output structure.


TODO: Resourcen in den parametern anpassen

TODO: Welche nextflow version soll es sein?

## Preparations

### Link to your volume
In this course, we are using the location `~/mgcourse/` as directory, where all the analysis is stored. We assume, that there is a larger storage volume mounted at `/vol/mgcourse`. If your prefered volume is at another location, change the following commands accordingly.

!!! question "Task 1"
    First, we will link the volume to `~/mgcourse`, so we can use that link for the rest of this tutorial:
    ```BASH
    sudo ln -s /vol/mgcourse/ ~/mgcourse
    ```

### Database directory link
In addition, we need a symlink for the directory where the databases should be stored.

!!! question "Task 2"
    Create a `databases` directory, the `/vol/scratch` directory and create a link to the database directory:
    ```BASH
    sudo mkdir ~/mgcourse/databases
    sudo mkdir -p /vol/scratch
    sudo ln -s ~/mgcourse/databases/ /vol/scratch/databases
    ```

### Download GTDB:
TODO: Anpassen fuer toolkit, oder einfach downloaden lassen beim toolkit aufruf? in der kaffeepause z.B.

First, we need to download the GTDB database files. The database is pretty
big (64 Gb), so even downloading it from our local copy in the de.NBI Cloud
will take a couple of minutes. After downloading, we need to extract the
tar archive (please be patient ;)::

  cd /mnt
  wget -qO- https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/denbi-mg-course/gtdbtk_v2_data.tar.gz | tar xvz
   
Now we need to set an environment variable that stores the path to
the database and put GTDBtk in our path::

  export GTDBTK_DATA_PATH=/mnt/release207_v2

  
We need to download the mash database, since it takes some time to create::

  cd /mnt/release207_v2
  wget https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/mg_databases/mash.msh

