import os,sys
sas_ids = sys.argv[1]
fo = open('rclone_command','w')
for sas_id in sas_ids.split(','):
    fo.write ('rm cal_solutions.h5.tar\n')
    command =  'timeout 900s rclone --config=/home/njj/macaroons/maca_sksp_tape_spiderlinc.conf copy maca_sksp_tape_spiderlinc:/L%s/cal_solutions.h5.tar ./'%sas_id
    fo.write (command+'\n')
    fo.write ('tar xvf cal_solutions.h5.tar\n')
    fo.write ('mv cal_solutions.h5 L%s_solutions.h5\n'%sas_id)

fo.close()
