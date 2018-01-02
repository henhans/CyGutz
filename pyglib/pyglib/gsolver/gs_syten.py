import sys, h5py, subprocess, os, numpy, glob

def write_text_2d(fname, head, dm):
    '''Write text file of 2d arrays.
    '''
    with open(fname, 'w') as f:
        f.write(head)
        for i, dm1 in enumerate(dm):
            for j, dm2 in enumerate(dm1):
                if abs(dm2) > 1.e-12:
                    f.write('{:3d}{:3d}{:24.16f}{:24.16f}\n'.format( \
                            i, j, dm2.real, dm2.imag))


def write_text_4d(fname, head, dm):
    '''Write text file of 4d arrays.
    '''
    with open(fname, 'w') as f:
        f.write(head)
        for i, dm1 in enumerate(dm):
            for j, dm2 in enumerate(dm1):
                for k, dm3 in enumerate(dm2):
                    for l, dm4 in enumerate(dm3):
                        if abs(dm4) > 1.e-12:
                            f.write('{:3d}{:3d}{:3d}{:3d}{:24.16f}{:24.16f}\n' \
                                    .format(i, j, k, l, dm4.real, dm4.imag))


def h5gen_text_embed_hamil(imp):
    '''Generate text files of the embedding Hamiltonian.
    '''
    with h5py.File('EMBED_HAMIL_{}.h5'.format(imp), 'r') as f:
        daalpha = f['/D'][...].T
        lambdac = f['/LAMBDA'][...].T
        h1e     = f['/H1E'][...].T
        v2e     = f['/V2E'][...].T

    # h1e file
    head = '# list of h1e_{alpha, beta} of significance \n' + \
            '# alpha   beta            h1e.real            h1e.imag \n'
    write_text_2d('h1e_{}.dat'.format(imp), head, h1e)

    # d file
    head = '# list of d_{a, alpha} with significance \n' + \
            '#     a  alpha              d.real              d.imag \n'
    write_text_2d('d_aalpha_{}.dat'.format(imp), head, daalpha)

    # lambdac file
    head = '# list of lambdac_{a, b} of significance \n' + \
            '#     a      b        lambdac.real        lambdac.imag \n'
    write_text_2d('lambdac_{}.dat'.format(imp), head, lambdac)

    # v2e file
    head = '# Note H_loc = \sum_{a,b} {h1e_a,b c_a^\dagger c_b} \n' + \
            '#            + 1/2 \sum_{a, b, c, d} v2e_{a,b,c,d}  \n' + \
            '#            * c_a^\dagger c_c^\dagger c_d c_b \n' + \
            '# list of v2e_{alpha, beta, gamma, delta} of significance \n' + \
            '# alpha   beta  gamma  delta              v2e.real' + \
            '              v2e.imag\n'
    if not os.path.isfile('v2e_{}.dat'.format(imp)):
        write_text_4d('v2e_{}.dat'.format(imp), head, v2e)


def gen_file_lat(imp, flat, norbs, reorder=None, threads_super=1):
    '''Generate lat file.
    '''
#    command = ['/usr/bin/time', '-v', '-a', '-o', 'timing_{}'.format(imp),
#            'syten-mps-embedding-lanata', '-l', str(norbs), '--h1e',
#	    './h1e_{}.dat'.format(imp), '--d', './d_aalpha_{}.dat'.format(imp), '--lambdac',
#            './lambdac_{}.dat'.format(imp), '--v2e', './v2e_{}.dat'.format(imp),
#            '-o', flat]
#    if reorder is not None:
#        command.extend(['-r', reorder])
#    command.extend(['--threads-super',str(threads_super), '-q'])
#    print(' '.join(command))
#    subprocess.call(command)
    cmd = '/usr/bin/time -v -a -o'+' timing_{}'.format(imp) + \
            ' syten-mps-embedding-lanata -l '+ str(norbs) + ' --h1e ' + \
	    ' ./h1e_{}.dat'.format(imp) + ' --d ' + ' ./d_aalpha_{}.dat'.format(imp) + ' --lambdac '+ \
            ' ./lambdac_{}.dat'.format(imp) + ' --v2e ' + ' ./v2e_{}.dat'.format(imp) + \
            ' -o ' + flat
    if reorder is not None:
        cmd += ' -r '+ reorder                          
    cmd += ' --threads-super ' + str(threads_super) #+ ' -q '
    print "Executing: ", cmd
    out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()


def gen_file_rnd_state(imp, norbs, lat, threads_tensor=1):
    '''Generate random state for initialization.
    '''
#    command = ['/usr/bin/time', '-v', '-a', '-o', 'timing_{}'.format(imp),
#            'syten-random', '-s', str(norbs), '-l', lat, '-o',
#            'rnd_{}.state'.format(imp), '--threads-tensor', str(threads_tensor), '-q']
#    print(' '.join(command))
#    subprocess.call(command)
    cmd = '/usr/bin/time -v -a -o timing_{}'.format(imp) + \
            ' syten-random -s ' + str(norbs) + ' -l ' + lat + ' -o ' + \
            ' rnd_{}.state'.format(imp) + ' --threads-tensor ' + str(threads_tensor) #+ ' -q '
    print "Executing: ", cmd
    out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()


def run_syten_dmrg(imp, lat, inp_state, out_f, s_config, out_state=None,
        threads_tensor=1):
    '''Run dmrg calculation.
    '''
#    command = ['/usr/bin/time', '-v', '-a', '-o', 'timing_{}'.format(imp),
#            'syten-dmrg', '-l', lat, '-i', inp_state, '-o', out_f,
#            '-s', s_config]
#    if out_state is not None:
#        command.extend(['-f', out_state])
#    command.extend(['--threads-tensor', str(threads_tensor), '-q'])
#    print(' '.join(command))
#    subprocess.call(command)
    cmd = '/usr/bin/time -v -a -o ' + ' timing_{}'.format(imp) + \
            ' syten-dmrg -l "' + lat + '" -i ' + inp_state + ' -o ' + out_f + \
            ' -s "' + s_config + '"'
    if out_state is not None:
        cmd += ' -f ' + out_state
    cmd += ' --threads-tensor ' + str(threads_tensor) #+ ' -q '
    print "Executing: ", cmd
    out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()


def run_syten_dmrg_tunnel(imp, lat, inp_state, out_f, s_config, out_state=None,
        threads_tensor=1):
    '''Run dmrg calculation with tunneling term to prevent stucking in local minimum.
    '''
#    command = ['/usr/bin/time', '-v', '-a', '-o', 'timing_{}'.format(imp),
#            'syten-dmrg', '-l', lat, '-l', 'with-reordering_{}.lat:Htunnel'.format(imp), '-i', inp_state, '-o', out_f,
#            '-s', s_config]
#    if out_state is not None:
#        command.extend(['-f', out_state])
#    command.extend(['--threads-tensor', str(threads_tensor)])#, '-q'])
#    print(' '.join(command))
#    log = open('syten-dmrg.log','w')
#    subprocess.call(command, stdout=log, stderr=log)
    cmd = '/usr/bin/time -v -a -o ' + ' timing_{}'.format(imp) + \
            ' syten-dmrg -l "' + lat + '" -l ' + ' "with-reordering_{}.lat:Htunnel"'.format(imp) + \
	    ' -i ' + inp_state + ' -o ' + out_f + ' -s "' + s_config +'" -v 4'
    if out_state is not None:
        cmd += ' -f ' + out_state
    cmd += ' --threads-tensor ' + str(threads_tensor)#+ ' -q '
    log = open('syten-dmrg.log','w')
    print "Executing: ", cmd
    subprocess.Popen(cmd, shell=True, stdout=log, stderr=log).communicate()


def run_syten_expectation(imp, lat, state, f_expval, threads_tensor=1):
    '''Calculate expectation values.
    '''
#    command = ['/usr/bin/time', '-v', '-a', '-o', 'timing_{}'.format(imp),
#            'syten-expectation', '-c', '-a', state, '-l', lat, \
#            '--template-file', f_expval, '--threads-tensor', str(threads_tensor), '-q']
#    return subprocess.check_output(command)
    cmd = '/usr/bin/time -v -a -o ' + ' timing_{}'.format(imp) + \
            ' syten-expectation -c -a ' + state + ' -l ' + lat + \
            ' --template-file ' + f_expval + ' --threads-tensor ' + str(threads_tensor) #+ ' -q '
    print "Executing: ", cmd
    out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    return out

def run_syten_mutual_information(imp, state, f_reorder, threads_tensor=1):
#    command = ['/usr/bin/time', '-v', '-a', '-o', 'timing_{}'.format(imp),
#            'syten-mutualInformation', '-o', state, '-f', f_reorder,
#            '--threads-tensor', str(threads_tensor), '-q']
#    print(' '.join(command))
#    subprocess.call(command)
    cmd = '/usr/bin/time -v -a -o ' + ' timing_{}'.format(imp) + \
            ' syten-mutualInformation -o ' + state + ' -f ' + f_reorder + \
            ' --threads-tensor ' + str(threads_tensor) #+ ' -q '
    print "Executing: ", cmd
    out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()

def driver_dmrg(s_config= \
        '(t 1e-8 d 1e-8 m 50 expandBlocksize 10 x 5 save false)' + \
        '(m 100) '*8 + '(m 100 d 1e-8 save true)'):
    '''driver for the impurity dmrg calculations.
    '''
    if '-i' in sys.argv:
        imp = sys.argv[sys.argv.index('-i')+1]
    else:
        imp = 1

    if os.path.isfile('s_config_{}.cfg'.format(imp)):
        s_config = ''
        for line in open('s_config_{}.cfg'.format(imp), 'r').readlines():
            s_config += line.replace('\n', '')

    if os.path.isfile('timing_{}'.format(imp)):
        os.remove('timing_{}'.format(imp))

    if os.path.isfile('s_threads_tensor.in'):
        threads_tensor = numpy.loadtxt('s_threads_tensor.in',dtype=int)
    else:
        threads_tensor =  1
    if os.path.isfile('s_threads_super.in'):
        threads_super = numpy.loadtxt('s_threads_super.in',dtype=int)
    else:
        threads_super = 1

    with h5py.File('EMBED_HAMIL_{}.h5'.format(imp), 'r') as f:
        na2 = f['/na2'][0]
        lambdac = f['/LAMBDA'][...].T

    # dump to text file due to dmrg code
    h5gen_text_embed_hamil(imp)

    # remove the previous temporary file
#    for filename in glob.glob('with*reordering*')+glob.glob('rnd*')+glob.glob('reordering*'):
#	os.remove(filename)
#    
#    # reorder the lattice
#    f_reorder = 'reordering_{}'.format(imp)
#    # generate lat file
#    flat = 'without-reordering_{}.lat'.format(imp)
#    gen_file_lat(imp, flat, na2, threads_super=threads_super)
#
#    # generate random state for initialization.
#    gen_file_rnd_state(imp, na2, flat, threads_tensor=threads_tensor)
#
#    # initial stage 1/2 sweeping
#    flat = 'without-reordering_{}.lat:H1e Hd Hf Hv2e + + +'.format(imp)
#    inp_state = 'rnd_{}.state'.format(imp)
#    out_f = 'without-reordering-dmrg_{}'.format(imp)
#    s_config0 = '(t 1e-16 d 1e-14 m 50 x 5) (m 100)'
#    run_syten_dmrg(imp, flat, inp_state, out_f, s_config0, threads_tensor=threads_tensor)
#
#    # reorder
#    state = out_f + '_2_5.state'
#    run_syten_mutual_information(imp, state, f_reorder, threads_tensor=threads_tensor)
#
#    if os.path.isfile('with-reordering-gs_{}.state'.format(imp)):
#        inp_state = 'with-reordering-gs_{}.state'.format(imp)
#    else:
#        inp_state = 'rnd_{}.state'.format(imp)

    # keep all the previous tempoarary files
    # reorder the lattice
    f_reorder = 'reordering_{}'.format(imp)
    if not os.path.isfile(f_reorder): 
        # generate lat file
    	flat = 'without-reordering_{}.lat'.format(imp)
    	gen_file_lat(imp, flat, na2, threads_super=threads_super)

    	# generate random state for initialization.
    	gen_file_rnd_state(imp, na2, flat, threads_tensor=threads_tensor)

    	# initial stage 1/2 sweeping
    	flat = 'without-reordering_{}.lat:H1e Hd Hf Hv2e + + +'.format(imp)
    	inp_state = 'rnd_{}.state'.format(imp)
    	out_f = 'without-reordering-dmrg_{}'.format(imp)
    	s_config0 = '(t 1e-16 d 1e-14 m 50 x 5) (m 100)'
    	run_syten_dmrg(imp, flat, inp_state, out_f, s_config0, threads_tensor=threads_tensor)

    	# reorder
    	state = out_f + '_2_5.state'
    	run_syten_mutual_information(imp, state, f_reorder, threads_tensor=threads_tensor)
    else:
        if os.path.isfile('with-reordering-gs_{}.state'.format(imp)):
            inp_state = 'with-reordering-gs_{}.state'.format(imp)
        else:
            inp_state = 'rnd_{}.state'.format(imp)

    # generate lat file with reordering
    flat = 'with-reordering_{}.lat'.format(imp)
    gen_file_lat(imp, flat, na2, reorder=f_reorder, threads_super=threads_super)

    # dmrg after reordering
    flat = 'with-reordering_{}.lat:H1e Hd Hf Hv2e + + +'.format(imp)
    out_f = 'with-reordering-dmrg_{}'.format(imp)
    out_state = 'with-reordering-gs_{}.state'.format(imp)
    run_syten_dmrg_tunnel(imp, flat, inp_state, out_f, s_config,
            out_state=out_state, threads_tensor=threads_tensor)

    # get expectation value
    f_temp = 'exp_val_{}.template'.format(imp)
    na4 = na2*2
    if not os.path.isfile(f_temp):
        with open(f_temp, 'w') as f:
            for i in range(na4):
                for j in range(i, na4):
                    f.write('{{ CH_{} C_{} * }}\n'.format(i, j))
            f.write('{ H1e Hd Hf Hv2e + + + }')
    flat = 'with-reordering_{}.lat'.format(imp)
    res = run_syten_expectation(imp, flat, out_state, f_temp, threads_tensor=threads_tensor)
    res = res.split('\n')[:na4*(na4+1)/2+1]
    dm = numpy.zeros([na4, na4], dtype=numpy.complex)
    ij = 0
    for i in range(na4):
        for j in range(i, na4):
            res_ = map(float, res[ij].split())
            ij += 1
            dm[i, j] = res_[0] + res_[1]*1.j
            if i != j:
                dm[j, i] = numpy.conj(dm[i, j])
    res_ = map(float, res[-1].split())
    assert numpy.abs(res_[1]) < 1.e-6, ' syten: non-real impurity etot!'

    # convention: f_b f_a^\dagger = \delta_{a,b} - f_a^\dagger f_b
    etot = res_[0] - numpy.trace(lambdac)
    with h5py.File('EMBED_HAMIL_RES_{}.h5'.format(imp), 'w') as f:
        f['/emol'] = [etot.real]
        f['/DM'] = dm.T



if __name__ == '__main__':
    driver_dmrg()
