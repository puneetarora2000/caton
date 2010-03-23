*********************************
Performance validation
*********************************

Intracellular-extracellular dataset
===========================================

We test the program on the data in which a cell is recorded, *in vivo*, intracellularly and extracellularly. (The extracellular recording also contains other cells.) See (Henze et al., 2000) for a description of the methods. The dataset can be downloaded at http://crcns.org/data-sets/hc/hc-1.

We run the clustering script on all files in the ``d11221`` group. We also extract the intracellular spikes. We compare these: an extracellular spike is matched to the intracellular spike if the time difference is less than .5 ms. We then compare the intracellular unit to the extracellular unit that matches it the most times. The following report is generated:

.. literalinclude:: generated_data/d11221/report.txt

In most cases, the intracellular unit produced a weak signal and was grouped together with other low-amplitude units into a big cluster. The exceptions are d11221.002.dat, d1122105.dat, d1122107.dat, d1122108.dat, and d1122109.dat

If you want to verify these results (or, better, change various parameters and see if you can get better results), download the ``d11221`` dataset and extract the zip archive. Navigate to ``caton/doc`` and change the line ``DataDir = ...`` to point to the correct location. Then run the script ``test_d11221.py``. The report will be in ``caton/doc/generated_data/d11221/report.txt``. This directory also contains the clustering results (.clu files, etc.)


Performance on synthesized 32-channel data
===========================================

Next, we do the only truly verifiable test we can do on 32-channel data. We take 4 channels of data from d11221.002.dat, and stack 8 copies of it with a time offset. Then we tell the program that this data comes from a 32-site probe with two rows of channels (same as example probe in :ref:`probefiles`.)

.. literalinclude:: generated_data/fake32/report.txt

The type II error rate is higher, presumably because there are overlapping spikes that we get with the stacked 32-channel that didn't occur with the original 32-channel data. The type I error rate is lower on some recordings--this may have to do with thresholds being different (threshold is determined by a segment of several seconds at the start of the data.)


Performance on real 32-channel data
====================================

Next we run the program on data collected from a 32-channel probe, recorded from the entorhinal cortex of a sleeping rat. The recording also contains electrical stimulation. This data was collected by Mariano Belluscio from the Buzs |a| ki lab.

The program detected 1 million spikes and grouped them into 24 clusters. Clusters 18-24 are artifacts from stimulation. The following are the mean waveforms and autocorrelograms of the remaining clusters:

.. image:: images/klusters1.png
   :width: 15cm
   :align: center

.. image:: images/klusters2.png
   :width: 15cm
   :align: center

.. image:: images/klusters3.png
   :width: 15cm
   :align: center

.. image:: images/klusters4.png
   :width: 15cm
   :align: center

.. image:: images/klusters5.png
   :width: 15cm
   :align: center

.. image:: images/klusters6.png
   :width: 15cm
   :align: center

.. image:: images/klusters7.png
   :width: 15cm
   :align: center

.. image:: images/klusters8.png
   :width: 15cm
   :align: center

.. image:: images/klusters9.png
   :width: 15cm
   :align: center

.. image:: images/klusters10.png
   :width: 15cm
   :align: center

.. image:: images/klusters11.png
   :width: 15cm
   :align: center

.. image:: images/klusters12.png
   :width: 15cm
   :align: center

.. image:: images/klusters13.png
   :width: 15cm
   :align: center

.. image:: images/klusters14.png
   :width: 15cm
   :align: center

.. image:: images/klusters15.png
   :width: 15cm
   :align: center

.. image:: images/klusters16.png
   :width: 15cm
   :align: center


.. image:: images/klusters17.png
   :width: 15cm
   :align: center

Last, here is the cross-correlogram of the non-artifact clusters:

.. image:: images/ccg.png
   :width: 20cm
   :align: center

.. |a| unicode:: U+00E1
   :trim:
