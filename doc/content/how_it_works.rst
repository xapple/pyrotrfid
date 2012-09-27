============
How it works
============

Mapping to reference
--------------------
The mapping is done using the BWA-SW algorithm and works for single-end reads only. BWA-SW is a tool designed for long reads. It is based on Burrows-Wheeler Transform (BWT). It performs heuristic Smith-Waterman-like alignment to find high-scoring local hits.

Digital digestion
-----------------
dT-RFLP stands for digital terminal-restriction fragment length polymorphism. For digital digestion, the software uses the ``biopython`` restriction library which contains restriction enzymes information and associated-methods to manipulate them as objects.

Flowchart
---------

.. image:: /images/flowchart.png