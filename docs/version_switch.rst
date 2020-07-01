.. _version:

==========================================================
Switching from version 1 (PyPopStar) to version 2 (SPISEA)
==========================================================
Version 1 of this software package was called PyPopStar. Due to a naming conflict,
we have renamed the package of SPISEA for version 2 and beyond.
For reference, v2.0.0 has the same functionality as v1.0.1. If you already downloaded
v1.0.0 or v1.0.1, the steps required to upgrade to v2.0.0 are described below. 


Step 1: Get Updated Package from Github
-------------------------------------
To get the updated package from Github, you need to pull down the updated code respository.
In your PyPopStar directory, run the following command from the terminal::
  git pull

[ADD INSTRUCTIONS: git checkout main]



Step 2: Change PYTHONPATH To Point to New Directory Name
----------------------------------------------------------




Step 3: Change Environment Variable Name
------------------------------------------
Update the name of the environment variable ``POPSTAR_MODELS`` to
the new name ``SPISEA_MODELS``. If you are using bash, this line in your
.bash_profile would be::
  export SPISEA_MODELS=/<path_to_models_directory>


Step 4: Test the Package to Make Sure It Is Working
---------------------------------------------------
To make sure everything is working, go through the instructions on
:ref:`test-setup`. 



Step 5: Change Import Statements in Subsequent Code
---------------------------------------------------
In all code that calls SPISEA functions (e.g. synthetic.py), you will need to
change the import statements from::
  from popstar import synthetic
  
to::
  from spisea import synthetic





