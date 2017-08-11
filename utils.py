import os
import sys
import time


def format_time(fp, header, str):
    fp.write('==' + time.strftime("%b %d %Y %H:%M:%S", time.localtime()) + '== [' + header + '] ' + str)


def exec_cmd(fp, header, cmd):
    format_time(fp, header, '[CMD] '+ cmd + '\n')
    os.system(cmd)
