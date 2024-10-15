## Database input configuration

Whenever a database field can be specified as part of the tool configuration (such as in gtdb or checkm), you are able to provide different methods to
fetch the database. In all settings, please make sure that the file has the same ending (e.g. .zip, .tar.gz) as specified in the corresponding tool section.
In addition, as database names are used to name results with which they were created, said database names should contain the respective database number or date of creation.
With this every result can be linked to one exact database version to clarify results. 
Except for the `extractedDBPath` parameter, all other input types (https, s3,...) will download the database to the folder specified in the `database` parameter.

### Extracted Database Path

If you have already downloaded and extracted the database, you can specify the path using the `extractedDBPath` parameter.
This setting is available in standard and slurm mode. In slurm mode the path can point to a db on the worker node.

Example:
```
database:
  extractedDBPath: /vol/spool/gtdb/release202
```

### HTTPS Download

The toolkit is able to download and extract the database, as long as the file ending equals the one specified in the corresponding tool section (.zip, tar.gz, tar.zst)
This setting is available in standard and slurm mode. 


Example:
```
database:
  download:
    source: 'https://openstack.cebitec.uni-bielefeld.de:8080/databases/gtdb.tar.gz'
    md5sum: 77180f6a02769e7eec6b8c22d3614d2e 
```

### Local File Path

This setting allows you to reuse an already downloaded database. 

Example:
```
database:
  download:
    source: '/vol/spool/gtdb.tar.gz'
    md5sum: 77180f6a02769e7eec6b8c22d3614d2e 
```

### S3 Download

You can specify an S3 link and configure the S3 call via the `s5cmd.params parameter.
The `s5cmd.params` parameter allows you to set any setting available of the [s5cmd](https://github.com/peak/s5cmd) commandline tool. 
If you need credentials to access your databases, you can set them via the Nextflow secrets mechanism. The correct key name for the access and secret key can be found in the corresponding database section.

In the following example the compressed file will be downloaded and extracted.

Example for publicly available compressed database:
```
database:
  download:
    source: 's3://databases/gtdb.tar.gz'
    md5sum: 77180f6a02769e7eec6b8c22d3614d2e 
    s5cmd:
      params: '--retry-count 30 --no-sign-request --no-verify-ssl --endpoint-url https://openstack.cebitec.uni-bielefeld.de:8080'
```

If your database is already extracted and available via S3, you can specify the S3 link using a wildcard as in the next example.

```
database:
  download:
    source: 's3://databases/gtdb/*'
    md5sum: 77180f6a02769e7eec6b8c22d3614d2e 
    s5cmd:
      params: '--retry-count 30 --no-verify-ssl --endpoint-url https://openstack.cebitec.uni-bielefeld.de:8080'
```

### Updating Database MD5SUMs 

The md5sum is computed over all md5sums of all files of the extracted database.
If you need to update the md5sum because you updated your database you have to download the database 
and run the following command

```
find /path/to/db -type f -exec md5sum {} \; | sort | cut -d ' ' -f 1 | md5sum | cut -d ' ' -f 1
```

### Database Download strategy

The toolkit allows to download databases on multiple nodes and tries to synchronize the download process between
multiple jobs on a node. However not all possible combinations of profiles and download types are reasonable.

| PROFILE     | Download to Shared NFS | Download to worker database directory | Reuse extracted directory |
| :---------- | :--------------- | :------------------------------------ | :------------------------ |
| STANDARD    | :material-check: |  :material-close:                     | :material-check:          |
| SLURM       | :material-close: | :material-check:                  | :material-check:  On a disk local to the worker and nfs directory |


