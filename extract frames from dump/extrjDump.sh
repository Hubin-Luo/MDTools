#!/bin/bash

#extract one dump frame. Usage: scriptname dumpfile iframe
#frame starts from 0

awk 'BEGIN{iframe=0}
     /TIME/{
        iframe+=1
        if(iframe=='$2'+1){
           print
           while(!/ITEM: ATOMS/){
              if(/NUMBER/){
                 getline; print; natom=$1
              }
              else{
                 getline; print
              }
           }
           if(/ITEM: ATOMS/){
              for(i=1;i<=natom;i++){
                 getline; print
              }
           }
        }
     }' $1 > extrj-$2.dump
