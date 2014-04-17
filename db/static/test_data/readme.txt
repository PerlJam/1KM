The test_data files can be used to initialize the database with a set of test data.

IMPORTANT: bootstrap files must be loaded before this data can be loaded.  see:
- reports/static/api_init/api_init_actions.csv 
- db/static/api_init/api_init_actions.csv 

- api_init_actions.csv is a recipe file used by the "db_init" script, telling
it what order to load each of the api_init files in the bootstrapping process. 
It is run like this:
run the server in a localhost port:
(virtualenv)$./manage.py runserver 55001
run the bootloader script:
(virtualenv)$ ./manage.py db_init  --inputDirectory=./db/static/test_data/ \
  -f ./db/static/test_data/api_init_actions.csv \
  -u http://localhost:55001/db/api/v1 -a <user:password>
