Reference
=========

.. currentmodule:: colorspacious

Colorspaces
-----------

.. ipython::

   import colorspacious
   with open("_static/colorspacious.dot", "w") as f:
       colorspacious.conversion.GRAPH.dump_dot(f)
   import subprocess
   subprocess.call(["dot", "-Tsvg", "_static/colorspacious.dot",
                    "-o", "_static/colorspacious.svg"])

.. image:: /_static/colorspacious.*

XX TODO

Functions and objects
---------------------

XX TODO
