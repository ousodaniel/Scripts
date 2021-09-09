## Bash practise

addgroup computestream

mkdir -p /exam/computestream

chgrp :computestream /exam/computestream

adduser 

#sudoer user privs: who where=(as_whom) what


adduser candidate

passwd candidate; cert456

visudo /etc/sudoers; candidate ALL=(ALL:ALL) NOPASSWORD:ALL

touch /etc/skel/NEWS

addgroup students

adduser  --home /home/school/harry/ --in-group students harry; passwd harry; magic; echo PATH=/home/school/harry/binaries:$PATH > \
/home/school/harry/.profile

adduser --home /sysadmin/ --shell /usr/bin/zsh; sysadmin ALL=(ALL:ALL) NOPASSWORD:ALL

User_Alias LAST=computestream,harry,sysadmin; LAST ALL=(ALL:ALL) /usr/bin/last

paaswd projectadmin; _onetime43_ | usermod --password _onetime43_ --home /home/projectadmin projectadmin

usermod --shell /usr/bin/sh devel

cat /etc/services | 2605 | cut -f1 > /home/student/port-2605.txt


cat /etc/services | grep imap | sed 's/\t\t/ /g' | grep -o "[0-9]*/" | tr -d / > /home/student/imap-ports.txt

mount /dev/xvdf2 /mnt/backup/

tar -jxvf /mnt/backup/backup-primary.tar.bz2 -C /opt/

/dev/xvdi1 none swap sw,noauto 0 0

/staging /mnt/point ext4 ro 0 0


ls -1  | xargs printf "$PWD/%s\n" > test.txt 

unzip /opt/SAMPLE001.zip -d /opt/SAMPLE001

tar /opt/SAMPLE001.tar /opt/SAMPLE001

bzip2 -z SAMPLE001.tar.bz /opt/SAMPLE001.tar