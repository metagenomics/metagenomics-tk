



Todos:
* Introduce the Toolkti
* Introduce the data to process
* Download gtdb!
* Set a working directory to a volume like ~/mgc_workdir/


## Download GTDB:
TODO: Anpassen fuer toolkit, brauchen wir mash?

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

