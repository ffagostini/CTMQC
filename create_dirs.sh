#!/bin/bash

dir=output
if [ ! -d $dir ] ; then
   mkdir $dir
else
  rm -rf $dir/*
fi

dir=output/histo
if [ ! -d $dir ] ; then
   mkdir $dir
else
  rm -rf $dir/*
fi

dir=output/coeff
if [ ! -d $dir ] ; then
   mkdir $dir
else
  rm -rf $dir/*
fi

dir=output/trajectories
if [ ! -d $dir ] ; then
   mkdir $dir
else
  rm -rf $dir/*
fi

dir=output/density
if [ ! -d $dir ] ; then
mkdir $dir
else
rm -rf $dir/*
fi




