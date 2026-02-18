=======================
Scipion AMBER plugin
=======================

In order to use this plug-in, you need to have `Scipion3 <https://scipion-em.github.io/docs/docs/scipion-modes/how-to-install.html>`_ installed.

Full documentation to this plugin can be found in the `official documentation page <https://scipion-chem.github.io/docs/plugins/amber/index.html>`_.

Then you have to follow the following steps:

**Clone this repository:**

.. code-block::

    git clone https://github.com/scipion-chem/scipion-chem-amber

**Install the plugin in devel mode**

.. code-block::

    scipion3 installp -p path_to_scipion-chem-amber --devel 

OR

**Install the plugin in user mode**

.. code-block::

    scipion3 installp -p path_to_scipion-chem-amber

**Binary files** 

AmberTools is automatically installed but it can only run simulations on CPU. 
To enable the GPU acceleration, non-comercial **Amber24 suite** is required.
You must either install pmemd in the EM_ROOT (typically: SCIPION_HOME/software/em/amber-202*-*)
or define the path of pmemd24 in scipion.conf as PMEMD_HOME.
For official installation instructions, visit: <https://ambermd.org/GetAmber.php#amber> 

