import math
import time
import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.mplot3d import Axes3D
import warnings
warnings.filterwarnings("ignore")

e  = 2.7182818
pi = 3.1415923

def disp(kcen, kwid, knum, p1, p2, x1, x2):
    
    ### Parameter Input
    
    nstep = 20000       # number of steps in x and k vectors
    nstd3 = 4           # number of standard deviations shown in 3D plot
    nstdk = 80          # number of standard deviations shown in 2D k plot
    nstdx = 4           # number of xmax lengths shown in 2D x plot
    
    # Note: for optimal uncertainty calculations, use largest nstep possible 
    # and smallest nstdk and nstdx values without cutting off non-zero wave.
    
    ###

    k   = []
    k3  = []
    xp  = []
    xn  = []
    phi = []
    ph3 = []
    tot = []
    knum = int(knum)
        
    for i in range(nstep+1):
        
        stepk = nstdk*kwid/nstep 
        stepx = nstdx*x2/nstep
        
        new = (kcen-(nstdk/2)*kwid)+i*stepk
        k.append(new)
        
        xval = i*stepx
        xp.append(xval)
        xn.append(-1*xval)
        
        phi.append( (1/(kwid*(2*pi)**.5))* (e**(-(new-kcen)**2/(2*kwid**2) )) )     
        tot.append(0.)
        
        kmini = (nstep/2)-(nstd3*nstep/nstdk)
        kmaxi = (nstep/2)+(nstd3*nstep/nstdk)
        if i >= kmini and i <= kmaxi:
            k3.append(k[i])
            ph3.append(phi[i])
    
    plt.figure()
    gs = gridspec.GridSpec(2,2)
    
    ax3 = plt.subplot(gs[:, 1], projection='3d')
    ax3.plot(k3, ph3, 0, zdir='y', linewidth='3.0')

    for i in range(knum):
        frq = (kcen-nstd3*kwid)+i*2*nstd3*kwid/(knum-1)
        amp = phi[int(kmini+((kmaxi-kmini)/(knum-1)*i))]
        y  = []
        y3 = []
        x3 = []
        
        for j in range(nstep+1):
            y.append(amp*math.cos(frq*xp[j]))
            tot[j] += y[j]
            if j <= (nstep/nstdx)+1:
                y3.append(y[j])
                x3.append(xp[j])
        ax3.plot(x3, y3, frq, zdir='x')
    
    zlim = 1/(2.1*kwid)
    ax3.set_xlim3d(kcen-nstd3*kwid, kcen+nstd3*kwid)
    ax3.set_ylim3d(-.1, x2)
    ax3.set_zlim3d(-zlim, zlim)
    ax3.set_zticks([0])
    
    ax3.set_xlabel("momentum space (k)")
    ax3.set_ylabel("position space (x)")
    ax3.set_zlabel("Relative Amplitude")

    fig1 = plt.subplot(gs[0, 0])
    plt.plot(k, phi)
    plt.ylabel("Relative Amplitude")
    fig1.set_yticks([0])
    plt.xlim((p1, p2))
    plt.title('momentum space (k)')
    
    fig2 = plt.subplot(gs[1, 0])
    plt.plot(xp, tot, color='blue')
    plt.plot(xn, tot, color='blue')
    plt.ylabel("Relative Amplitude")
    fig2.set_yticks([0])
    plt.xlim((x1, x2))
    plt.title('position space (x)')
    
    intx = 0
    intk = 0
    for i in range(nstep+1):
        intx += 2*stepx*(tot[i]**2)
        intk +=   stepk*(phi[i]**2)

    sigx = 0
    sigk = 0
    for i in range(nstep+1):
        sigx += (2 * (xp[i]**2))*(tot[i]**2)*stepx/intx
        sigk += ((kcen-k[i])**2)*(phi[i]**2)*stepk/intk
        
    sigxk = round((sigx*sigk)**0.5, 3)
    sigx  = round(sigx**0.5, 3)
    sigk  = round(sigk**0.5, 3)
    
    fig1.annotate(r'$\sigma_k=%s$'%(sigk), xy=(p1,0), xytext=(.7,.88), textcoords='axes fraction')
    fig2.annotate(r'$\sigma_x=%s$'%(sigx), xy=(0,0), xytext=(.7,.88), textcoords='axes fraction')
    fig2.annotate(r'$\sigma_x \sigma_p =%s \hbar$'%(sigxk), xy=(0,0), xytext=(.5, -.23), textcoords='axes fraction', horizontalalignment='center', fontsize=16)
    
    #fig1.hline(0, xmin=kcen-sigk, xmax=kcen+sigk, color='red', linewidth='3.0')
    #fig2.hline(0, xmin=-sigx, xmax=sigx, color='red', linewidth='3.0')
    
    plt.show()
    plt.draw()
    
    stillopen = True
    while (stillopen):
        plt.pause(.1)
        if not plt.get_fignums():
            stillopen = False
    return

def reset(orig, new):
    try:
        return float(new)
    except ValueError:
        print("\nWoah! That's not a number!")
        time.sleep(1)
        return orig

def main():
    
    print("\n")
    print("    Hello there! Welcome to")
    print("   Fourier Transformations &")
    print("   Understanding Uncertainty")
    print("\n")
    
    cen = 10
    wid = 1
    num = 11
    klo = 5
    khi = 15
    xlo = -5
    xhi = 5
    done = False
        
    while(not done):
        
        print("\nYour current transformation is:\n")
        print("Gaussian k-distribution centered at", cen, "with sigma", wid)
        print("showing", num, "component waves,", klo, "< k <", khi, "&", xlo, "< x <", xhi, "\n")        
        
        print("To view your current transformation, just type \"go\".")
        print("To adjust the Gaussian center, type \"cen=#\" for center at #.")
        print("To adjust the Gaussian sigma, type \"sig=#\" for a sigma of #.")
        print("To adjust the number of wave components shown, type \"num=#\".")
        print("To adjust k scale from a to b, type \"klo=a\" and \"khi=b\".")
        print("To adjust x scale from c to d, type \"xlo=c\" and \"xhi=d\".")
        print("Separate multiple commands with a comma.")
        print("If you wanna quit, just say \"bye\".")
        
        inp = input("\n>>> ")
        inps = inp.split(',')
        for j in range(len(inps)):
            inps[j] = inps[j].strip()
        
        for i in range(len(inps)):  
            
            pair = inps[i].split('=')
            for j in range(len(pair)):
                pair[j] = pair[j].strip()
            
            if pair[0] == "go":
                disp(cen, wid, num, klo, khi, xlo, xhi)
            
            elif pair[0] == "cen":
                cen = reset(cen, pair[1])
                
            elif pair[0] == "sig":
                wid = reset(wid, pair[1])
                
            elif pair[0] == "num":
                num = reset(num, pair[1])
                
            elif pair[0] == "klo":
                klo = reset(klo, pair[1])
                
            elif pair[0] == "khi":
                khi = reset(khi, pair[1])
                
            elif pair[0] == "xlo":
                xlo = reset(xlo, pair[1])
                
            elif pair[0] == "xhi":
                xhi = reset(xhi, pair[1])
                
            elif pair[0] == "bye":
                print("\nSee ya!\n")
                time.sleep(1)
                done = True
                
            else:
                print("\nYou fool! Entirely invalid...\n")
                time.sleep(1)

    return
main()