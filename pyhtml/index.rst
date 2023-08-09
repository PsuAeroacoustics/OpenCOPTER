.. OpenCOPTER documentation master file, created by
   sphinx-quickstart on Wed Jul 13 16:08:17 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to OpenCOPTER's Documentation!
===================================

Introduction
============

OpenCOPTER is a library for quickly simulating the aerodynamics of multirotor vehicles with enough resolution to compute the sound produced by the vehicle. This library is written in the `D Programming Language <https://dlang.org/>`_ but is usable from Python for easy integration into the rest of the Python scientific computing ecosystem. This is the implementation of the work described by Rau and Greenwood :cite:p:`rau2022dynamic`. This user manual aims to cover both the theory behind the code as well as guidance on usage of the library.

The original goal of this code was to compute the aeroacoustics of eVTOL vehicles which typically use rigid rotors. As such it contains little in the way of rotor and blade dynamics such as flapping. However, it does contain inputs for blade flapping and will integrate it into the aerodynamics appropriately.

Getting Started
===============

.. toctree::
   :caption: Getting Started

   gettingstarted

Using the Code
==============

.. toctree::
   :caption: Using the Code

   using

Theory
======

.. toctree::
   :caption: Theory

   theory

.. _API Documentation:

API Documentation
=================

.. toctree::
   :caption: API Documentation

   apidoc


.. only:: html

   Indices and tables
   ==================

   * :ref:`genindex`
   * :ref:`search`
