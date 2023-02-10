# !/bin/bash

if [ -f /mnt/hdd/public/share_resource/OS_WES/cfDNA/.bashrc ];
then
   Loc=`find "$(pwd)" -name "CNV_annotation.sh"`;
   echo "alias CNA_annotation='bash ${Loc}'" >> /mnt/hdd/public/share_resource/OS_WES/cfDNA/.bashrc;
   echo "Installation is already finished";
else
   echo "Your .bashrc does NOT exist";
fi