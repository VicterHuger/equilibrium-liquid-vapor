## Create Virtual Enviroment
- Reference: https://code.visualstudio.com/docs/python/environments#_creating-environments

## Activate Virtual Enviroment
- This should be done everytime before starting coding

### Windows

```powershell

# Activate virtual environment on windows

# From command prompt
.venv\bin\activate.bat

# From power shell
.venv\Scripts\Activate.ps1

```

### Linux

```bash
source directory-path/.venv/bin/activate
```

## Recommended Extensions (VSCODE)
- njpwerner.autodocstring
- ms-python.autopep8
- ms-python.pylint
- ms-python.python
- KevinRose.vsc-python-indent

## Install Dependencies
- To install all dependencies of the project:
```bash
pip install -r requirements.txt
```

- To add any new  dependency, just do:

```bash
pip install name_of_dependency
```

- After this update requirements.txt package:
```bash
pip freeze > requirements.txt
```


## Tortoise ORM

### Initialization:

- For creating migration directory in root directory, run the command.
```bash
 aerich init -t database.tortoise.tortoise_config.TORTOISE
```

- First migrattion (only needed once, first person in the project to do this)
```bash

aerich init-db

```

### Generate Migration:
- Generate a new migration script for the changes in your models. 
```bash
aerich migrate --name your_migration_name
```

### Apply Migration:
- Apply the migrations to update the database schema
```bash
aerich upgrade
```

### Revert Migrations:
- If needed, you can revert a migration:
```bash
aerich downgrade
```

## Enviroment Variables
```bash
#Format example
DATABASE_URL=postgres://user_name:password@localhost:5432/db-name

```