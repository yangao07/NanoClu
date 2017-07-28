import os
import sys
import time


def format_time(header, str):
    sys.stderr.write('==' + time.strftime("%b %d %Y %H:%M:%S", time.localtime()) + '== [' + header + '] ' + str)


def exec_cmd(header, cmd):
    format_time(header, '[CMD] '+ cmd + '\n')
    os.system(cmd)
