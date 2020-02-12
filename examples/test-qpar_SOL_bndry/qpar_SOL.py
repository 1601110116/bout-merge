
savefig = 1

parts = [
    # 'SOL',
    'bndry_c',
    ]

if 'SOL' in parts:
    cases = 'data'
    qpar_n = 'qpar_nonlocal_3d_20'
    qpar_l = 'qpar_landau_3d1_20'
    qparn = collect(qpar_n, path=cases, yguards=True)
    qparl = collect(qpar_l, path=cases, yguards=True)
    x1 = np.arange(66, 130)
    x2 = np.arange(66, 132)
    xcore = 15
    xsol = 25

    plt.figure('qpar_SOL')
    plt.clf()
    plt.plot(x2, qparl[xsol, np.r_[x2], 0], 'o-', label='sol-l', mfc="None")
    plt.plot(x1, qparl[xcore, np.r_[x1], 0], 's-', label='core-l', mfc="None")
    plt.plot(x2, qparn[xsol, np.r_[x2], 0], '+-', label='sol-n')
    plt.plot(x1, qparn[xcore, np.r_[x1], 0], 'x-', label='core-n')
    plt.legend()
    plt.axhline(0, color='r', ls='--', lw=1)

if 'bndry_c' in parts:
    cases = 'data'
    qn2d = collect('qpar_nonlocal_0', path=cases, yguards=True)
    qn = collect('qpar_nonlocal_3d', path=cases, yguards=True)
    qn0 = collect('qpar_nonlocal_3d_0', path=cases, yguards=True)
    qn1 = collect('qpar_nonlocal_3d_1', path=cases, yguards=True)
    qn2 = collect('qpar_nonlocal_3d_2', path=cases, yguards=True)
    qn3 = collect('qpar_nonlocal_3d_3', path=cases, yguards=True)
    qn4 = collect('qpar_nonlocal_3d_4', path=cases, yguards=True)
    ql0 = collect('qpar_landau_3d_0', path=cases, yguards=True)
    ql1 = collect('qpar_landau_3d_1', path=cases, yguards=True)
    ql2 = collect('qpar_landau_3d_2', path=cases, yguards=True)
    ql3 = collect('qpar_landau_3d_3', path=cases, yguards=True)
    ql4 = collect('qpar_landau_3d_4', path=cases, yguards=True)
    x1 = np.arange(66, 130)
    x2 = np.arange(66, 132)
    xcore = 15
    xsol = 25

    plt.figure('qpar_SBC')
    plt.clf()
    """
    plt.plot(x2, ql0[xsol, np.r_[x2], 0], 'o-', label='sol-l', mfc="None")
    plt.plot(x1, ql0[xcore, np.r_[x1], 0], 's-', label='core-l', mfc="None")
    plt.plot(x2, qn0[xsol, np.r_[x2], 0], '+-', label='sol-n')
    plt.plot(x1, qn0[xcore, np.r_[x1], 0], 'x-', label='core-n')
    """
    plt.plot(x2, ql0[xsol, np.r_[x2], 0], 'o-', label='w/o coll, 0', mfc="None")
    plt.plot(x2, ql1[xsol, np.r_[x2], 0], 's-', label='w/o coll, 1', mfc="None")
    plt.plot(x2, ql2[xsol, np.r_[x2], 0], '^-', label='w/o coll, 2', mfc="None")
    plt.plot(x2, qn0[xsol, np.r_[x2], 0], '-', label='w/ coll, 0', mfc="None")
    plt.plot(x2, qn1[xsol, np.r_[x2], 0], '-', label='w/ coll, 1', mfc="None")
    plt.plot(x2, qn2[xsol, np.r_[x2], 0], '-', label='w/ coll, 2', mfc="None")
    plt.legend(title='bndry_c flag') # fontsize=22)
    plt.axhline(0, color='r', ls='--', lw=1)
    plt.xlabel('Y/cm')
    plt.title('parallel heat flux(LF)')
    if savefig:
        bv.savefig('qpar_SBC_flag')

    plt.figure("qpar_SBC34")
    plt.plot(x2, ql0[xsol, np.r_[x2], 0], 'o-', label='w/o coll, 0', mfc="None")
    plt.plot(x2, ql1[xsol, np.r_[x2], 0], 's-', label='w/o coll, 1', mfc="None")
    plt.plot(x2, ql3[xsol, np.r_[x2], 0], '*-', label='w/o coll, 3', mfc="None")
    plt.plot(x2, ql4[xsol, np.r_[x2], 0], 'd-', label='w/o coll, 4', mfc="None")
    plt.plot(x2, qn0[xsol, np.r_[x2], 0], '-', label='w/ coll, 0', mfc="None")
    plt.plot(x2, qn1[xsol, np.r_[x2], 0], '-', label='w/ coll, 1', mfc="None")
    plt.plot(x2, qn3[xsol, np.r_[x2], 0], '-', label='w/ coll, 3', mfc="None")
    plt.plot(x2, qn4[xsol, np.r_[x2], 0], '-', label='w/ coll, 4', mfc="None")
    plt.legend(title='bndry_c flag') # fontsize=22)
    plt.axhline(0, color='r', ls='--', lw=1)
    plt.xlabel('Y/cm')
    plt.title('parallel heat flux(LF)')
    if savefig:
        bv.savefig('qpar_SBC_flag34')


plt.show()
