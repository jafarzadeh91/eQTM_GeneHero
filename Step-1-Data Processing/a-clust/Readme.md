Readme
================
Sina Jafarzadeh
2018-04-02

This directory belongs to R source code for a-clust algorithm. The two first files was written by us to use A-clust on our dataset using different appropriate settings and generate the reduced array of methylation probes.

clustering\_methylation\_probes\_different\_settings.R
======================================================
In this file, we run the a-clust algorithm for different settings of parameters. This file is considered the main file in this folder.

clustering\_methylation\_probes.R
=================================
In this file, we cluster methylation probes according the parameters set in the file above.

Acluster.R
==========
The main interface to a-clust algorithm. Files below are the internal functions utilized by this function and written by a-clust developers.

calc.dist.clust.point.R
=======================

calc.dist.clusters.R
====================

calc.dist.d.neighbor.R
======================

Dbp.merge.R
===========

first.R
=======

last.R
======

update.clust.indicator.R
========================

base.log
========
