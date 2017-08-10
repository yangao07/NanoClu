import threading as td
from os import listdir
from os import path

import utils


def poa_thread(poa_in_n, run_poa_sh, poa, fxtools, in_list, poa_mat):
    global POA_THREAD_I
    global POA_THREAD_LK

    while 1:
        POA_THREAD_LK.acquire()
        list_i = POA_THREAD_I
        POA_THREAD_I += 1
        POA_THREAD_LK.release()

        if poa_in_n <= list_i:
            break

        poa_in_list = in_list[list_i]

        clu_in = poa_in_list[0]
        base_name = poa_in_list[1]
        clu_out = poa_in_list[2]
        clu_out_tmp = poa_in_list[3]

        utils.exec_cmd('run_poa', 'bash %s %s %s %s %s %s %s %s'
                       % (run_poa_sh, poa, fxtools, clu_in, clu_out_tmp, poa_mat, base_name, clu_out))
    return


def run_poa(poa, thread_n, clu_folder, cons):
    global POA_THREAD_I
    global POA_THREAD_LK

    run_poa_sh = './run_poa.sh'
    poa_mat = '/u/home/y/yangao/software/poaV2/blosum80.mat'
    fxtools = '/u/home/y/yangao/bin/non_dot_fxtools'
    poa_in_list = []
    poa_in_n = 0
    num_list = listdir(clu_folder)

    threads = []
    POA_THREAD_LK = td.Lock()

    for num_dir in num_list:
        num_dir = clu_folder + num_dir + '/'
        fn_list = listdir(num_dir)
        for clu_in in fn_list:
            base_name = path.basename(clu_in)
            clu_in = num_dir + clu_in
            clu_out = clu_in + '.cons.fa'
            clu_out_tmp = clu_in + '.cons.tmp'

            poa_in_list.append((clu_in, base_name, clu_out, clu_out_tmp))
            poa_in_n += 1

    POA_THREAD_I = 0
    for i in range(thread_n):
        t = td.Thread(target=poa_thread, args=(poa_in_n, run_poa_sh, poa, fxtools, poa_in_list, poa_mat,))
        threads.append(t)
        t.start()
    for t in threads:
        t.join()

    utils.exec_cmd('run_poa', 'cat %s/*/*.cons.fa > %s' % (clu_folder, cons))
    utils.exec_cmd('run_poa', 'rm %s/*/*.cons.fa' % clu_folder)
    return
