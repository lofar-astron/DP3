 


Useful Commands
============================================

This page contains sets of command that can be used in combination with DP3 or as preparation on measurement sets before handing them over to DP3.


Creating Single Baseline MS
************************************

The following commands can be used to reduce a measurement set to a single baseline and stretch its visibilities over time.

The antenna names can be retrieved using:

.. code-block:: sh

  taql 'SELECT * FROM in.MS::ANTENNA'

You can then run:

.. code-block:: sh

  # Extract the baseline with antenna pair CS002HBA0 and RS305HBA.
  DPPP msin=in.MS msout=temp.MS steps=[filter] filter.baseline='CS002HBA0&RS305HBA' filter.remove=True

  # Optionally you can check if the size was reduced
  du -hs *.MS

  # Delay the TIME of the last visiblity so that there will be mising time slots
  taql 'UPDATE temp.MS SET TIME=(select TIME+100*INTERVAL from ::)[0] WHERE TIME IN (select gmax(TIME) from ::)'

  # THIS DOES NOT WORK YET
  taql 'UPDATE temp.MS SET FLAG=false'

  # Fill the measurement set with visibilites for the "missing" time slots
  DPPP msin=temp.MS msout=out.MS steps=[]

  # Remove the temporary file
  rm -rf ./temp.MS

Building documentation locally
************************************

If you want to generate the documentation from a cloned repository you will first need the following dependencies.
Alternatively Ninja can be used to build the docs as is shown in the gitlab-ci.yml file.

.. code-block:: sh

  pip3 install jsonschema2rst sphinx sphinx-rtd-theme 

Before building the docs you will need to convert the yml schemes do rst files:

.. code-block:: sh

  jsonschema2rst ../DP3/docs/schemas/ ../DP3/docs/steps/
  sphinx-build -b html ../DP3/docs docs

If you have built the documentation using ssh you can mount it on your local machine to open the html files:

.. code-block:: sh

  sudo sshfs -o allow_other,default_permissions,IdentityFile=~/.ssh/id_rsa user@host:/path/to/dp3/build/docs/ /local/machine/path
