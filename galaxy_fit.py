#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Práctica descomposición fotométrica - Física Extragaláctica

"""
import subprocess
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
import os
from photutils.isophote import EllipseGeometry
from photutils.isophote import Ellipse
from photutils.aperture import EllipticalAperture
import shutil
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.special import gamma
from astropy.cosmology import FlatLambdaCDM


############################## Archivos ###################################
galaxia = "NGC6132_i.fits" # Imagen original
mask = "NGC6132_mask2D_new.fits"
psf = "mofeta.fits"

config_exp = "config_exponential.dat"
config_model = "config_model.dat" # Archivo de configuración para ajustar modelos de exponencial + sersic
config_psf = "config_makeimage_moffat_psf.dat" # Archivo de configuración para la moffattiana


###################### Parámetros de la observación ######################
sky = 177.18507 # Cuentas de cielo
gain = 4.8850000 # Ganancia
ron = 6.2500000 # Readout noise
zcal = -23.711000 # Constante de calibración a magnitudes
pixel_size = 0.396 # arcsec/pixel

redshift = 0.01657 # datos de CDS
cosmo = FlatLambdaCDM(H0=67.4, Om0=0.3)
scale = (1./cosmo.arcsec_per_kpc_proper(0.01657)).value

# Secondary X-axis
def arcsec_kpc(x):
    return x * scale  # Example conversion: double the x values

def kpc_arcsec(x):
    return x / scale  # Example conversion: double the x values

# Secondary X-axis
def pix_kpc(x):
    return x *pixel_size * scale  # Example conversion: double the x values

def kpc_pix(x):
    return x / scale / pixel_size  # Example conversion: double the x values
############################ LaTeX rendering ##############################
plt.rc('text', usetex=True)
plt.rc('font', family='sans-serif')  # Use a serif font for LaTeX rendering
plt.rc('font', size=16)  # Adjust size to your preference
# Define the LaTeX preamble with siunitx
plt.rcParams['text.latex.preamble'] = r'''
            \usepackage{siunitx}
            \sisetup{
              detect-family,
              separate-uncertainty=true,
              output-decimal-marker={.},
              exponent-product=\cdot,
              inter-unit-product=\cdot,
            }
            \DeclareSIUnit{\cts}{cts}
            \DeclareSIUnit{\dyn}{dyn}
            \DeclareSIUnit{\mag}{mag}
            \DeclareSIUnit{\arcsec}{arcsec}
            \DeclareSIUnit{\parsec}{pc}
            \usepackage{sansmath}  % Allows sans-serif in math mode
            \sansmath
            '''
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": "Computer Modern Serif",
})

def ver(imagen):
    """
    Función para visualizar una imagen con un mapa de color

    Parameters
    ----------
    imagen : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    


if __name__ == "__main__":
    
    # Carpeta donde vamos a guardar plots generales:
    if os.path.exists("plots_generales"):
        shutil.rmtree("plots_generales")  
    os.makedirs("plots_generales") 
    
    ######################################################################
    # Imagen original con cuentas en cada píxel:
    ######################################################################
    hdu_list_galaxy = fits.open(galaxia)
    original = hdu_list_galaxy[0].data
    hdu_list_galaxy.close()
    # Imagen original sin cuentas negativas:
    #original_abs = abs(original)
    minimo_pos = np.min(original[original>0])
    original_abs=original
    original_abs[original_abs<0]=0.0


    ######################################################################
    # Imagen calibrada en mag/arcsec^2:
    ######################################################################
    calibrada = -2.5*np.log10(original/pixel_size**2)-zcal
    calibrada_abs = -2.5*np.log10(original_abs/pixel_size**2)-zcal
    calibrada_abs = np.nan_to_num(calibrada_abs, nan=30, posinf=30, neginf=30)
    calibrada_abs = np.clip(calibrada_abs, 0, None)
    #print(calibrada_abs)
    ######################################################################
    # Plot de la imagen calibrada
    ######################################################################
    x_extent = np.arange(0, calibrada_abs.shape[1]) * pixel_size
    y_extent = np.arange(0, calibrada_abs.shape[0]) * pixel_size
    fig1, ax1 = plt.subplots(figsize=(14, 12))
    cax = ax1.imshow(calibrada_abs, cmap='inferno', origin='lower', aspect="equal", extent=[x_extent.min(), x_extent.max(), y_extent.min(), y_extent.max()], vmax=23)
    cbar = plt.colorbar(cax, ax=ax1, pad=0.13)
    cbar.set_label(r"$\mu$ [$\unit{\mag\per\arcsec\squared}$]", fontsize=30)
    cbar.ax.tick_params(labelsize=24)
    ax1.set_xlabel(r"$X$ [$\unit{\arcsec}$] ", fontsize=30)
    ax1.set_ylabel(r"$Y$ [$\unit{\arcsec}$]", fontsize=30)
    ax1.tick_params(axis='both', which='major', labelsize=24)
    
    kpcax = ax1.secondary_xaxis('top', functions=(arcsec_kpc, kpc_arcsec))
    kpcax.set_xlabel(r"$X$ [$\unit{\kilo\parsec}$] ", fontsize=30, labelpad=14)
    kpcax.tick_params(axis='both', which='major', labelsize=24)
    
    kpcay = ax1.secondary_yaxis('right', functions=(arcsec_kpc, kpc_arcsec))
    kpcay.set_ylabel(r"$Y$ [$\unit{\kilo\parsec}$]", fontsize=30, labelpad=14)
    kpcay.tick_params(axis='both', which='major', labelsize=24)
    fig1.savefig(os.path.join("plots_generales","imagen_calibrada.pdf"), format='pdf', bbox_inches='tight')
    plt.close(fig1)
    
    # fig2, ax2 = plt.subplots(figsize=(10, 10))
    # cax = ax2.imshow(calibrada_abs, cmap='inferno', origin='lower', aspect="equal")
    # cbar = plt.colorbar(cax, ax=ax2)
    # cbar.set_label(r"$\mu$ / $\unit{\mag\per\arcsec\squared}$", fontsize=30)
    # cbar.ax.tick_params(labelsize=24)
    # ax2.set_xlabel(r"$X$ / pixel ", fontsize=30)
    # ax2.set_ylabel(r"$Y$ / pixel", fontsize=30)
    # ax2.tick_params(axis='both', which='major', labelsize=24)
    # plt.show()
    ######################################################################
    # Generar la PSF
    ######################################################################
    comando = (
            f'cd {os.getcwd()} && '
            f'makeimage {config_psf} --output mofeta.fits'
            )
    result=subprocess.run(comando, shell=True, capture_output=True)
    # print(result.stderr)
    
    ######################################################################
    # Plot de la PSF (Moffattiana)
    ######################################################################
    hdu_list_PSF = fits.open(psf)
    psf_image = hdu_list_PSF[0].data
    fig3, ax3 = plt.subplots(figsize=(16, 12))
    cax = ax3.imshow(psf_image, cmap='inferno', origin='lower', aspect="equal", norm=PowerNorm(gamma=0.3))
    cbar = plt.colorbar(cax, ax=ax3)
    cbar.set_label(r"Intensidad Normalizada", fontsize=30)
    cbar.ax.tick_params(labelsize=24)
    ax3.set_xlabel(r"$X$ [px]", fontsize=30)
    ax3.set_ylabel(r"$Y$ [pc]", fontsize=30)
    ax3.tick_params(axis='both', which='major', labelsize=24)
    fig3.savefig(os.path.join("plots_generales","psf_mofeta.pdf"), format='pdf', bbox_inches='tight')
    plt.close(fig3)
    
    ######################################################################
    # Perfil de brillo
    ######################################################################
    # Máscara:
    hdu_list_mask = fits.open(mask)
    mask_image = hdu_list_mask[0].data
    hdu_list_mask.close()
    # Imagen original enmascarada:
    masked_original = original - original*mask_image
    # Imagen calibrada enmascarada:
    masked_calibrada = calibrada_abs.copy()
    masked_calibrada[mask_image==1]=0

    
    # Fit de isofotas sobre la original enmascarada (para no tener en cuenta estrellas):
    geometry = EllipseGeometry(x0=187, y0=265, sma=50, eps=0.7, pa=-75.0 * np.pi / 180.0) #Elipse inicial
    ellipse = Ellipse(masked_original, geometry)
    isolist = ellipse.fit_image(maxgerr=0.4, maxsma=120, step=0.1, linear=False)
    resultados = isolist.to_table()
    
    # Plot de las isofotas sobre la imagen calibrada (sobre la original queda más feo)
    fig4, ax4 = plt.subplots(figsize=(10, 10))
    cax = ax4.imshow(masked_calibrada , cmap='inferno', origin='lower', aspect="equal",
                     vmax=23, vmin=11)
    cbar = plt.colorbar(cax, ax=ax4, pad=0.16)
    cbar.set_label(r"$\mu$ [$\unit{\mag\per\arcsec\squared}$]", fontsize=30)
    cbar.ax.tick_params(labelsize=24)
    ax4.set_xlabel(r"$X$ [px] ", fontsize=30)
    ax4.set_ylabel(r"$Y$ [px]", fontsize=30)
    ax4.set_xlim(120,250)
    ax4.set_ylim(135,395)
    ax4.tick_params(axis='both', which='major', labelsize=24)
    
    kpcax = ax4.secondary_xaxis('top', functions=(pix_kpc, kpc_pix))
    kpcax.set_xlabel(r"$X$ [$\unit{\kilo\parsec}$] ", fontsize=30, labelpad=14)
    kpcax.tick_params(axis='both', which='major', labelsize=24)
    
    kpcay = ax4.secondary_yaxis('right', functions=(pix_kpc, kpc_pix))
    kpcay.set_ylabel(r"$Y$ [$\unit{\kilo\parsec}$]", fontsize=30, labelpad=14)
    kpcay.tick_params(axis='both', which='major', labelsize=24)
    
    for iso in isolist:
        if iso.sma > 0:
            aper = EllipticalAperture((iso.x0, iso.y0), iso.sma,
                              iso.sma * (1 - iso.eps),iso.pa)
            aper.plot(color='white')
    fig4.savefig(os.path.join("plots_generales","isofotas.pdf"), format='pdf', bbox_inches='tight')
    plt.close(fig4)
    
    # Creamos la curva de brillo superficial (pixeles, cuentas):
    smas = resultados["sma"]
    intensities = resultados["intens"]
    int_errors = resultados["intens_err"]
    
    fig5, ax5 = plt.subplots(figsize=(10, 10))
    ax5.plot(smas, intensities, "o", color="purple", markersize=9, alpha=0.7)
    ax5.set_xlabel(r"semieje mayor [pixel] ", fontsize=30)
    ax5.set_ylabel(r"cuentas", fontsize=30)
    ax5.tick_params(axis='both', which='major', labelsize=24)
    fig5.savefig(os.path.join("plots_generales","brillo_sin_calibrar.pdf"), format='pdf', bbox_inches='tight')
    plt.close(fig5)
    
    
    # Curva de brillo superficial calibrada (arcsecs, magnitudes):
    smas_cal = smas*pixel_size
    mu = -2.5*np.log10(intensities/pixel_size**2)-zcal
    # Fit lineal a la región exponencial:
    exponential_zone = smas_cal > 17
    slope, intercept = np.polyfit(smas_cal[exponential_zone], mu[exponential_zone], 1)
    mu_fit = np.polyval([slope, intercept], smas_cal)
    h=(2.5)/(np.log(10)*slope*pixel_size)
    I_0=10**((2.5*np.log10(pixel_size**2)-intercept-zcal)/(2.5))
    print("Ajuste en la región exponencial:")
    print(f"I_0 = {I_0}")
    print(f"h = {h}")
    fig6, ax6 = plt.subplots(figsize=(10, 10))
    ax6.plot(smas_cal, mu, "o", color="purple", markersize=9, alpha=0.7)
    ax6.plot(smas_cal, mu_fit, "--", color="darkgreen", linewidth=3, alpha=0.7)
    ax6.set_xlabel(r"semieje mayor [$\unit{\arcsec}$] ", fontsize=30)
    ax6.set_ylabel(r"$\mu$ [$\unit{\mag\per\arcsec\squared}$]", fontsize=30)
    ax6.tick_params(axis='both', which='major', labelsize=24)
    ax6.invert_yaxis()
    
    kpcax = ax6.secondary_xaxis('top', functions=(arcsec_kpc, kpc_arcsec))
    kpcax.set_xlabel(r"semieje mayor [$\unit{\kilo\parsec}$] ", fontsize=30, labelpad=14)
    kpcax.tick_params(axis='both', which='major', labelsize=24)

    
    fig6.savefig(os.path.join("plots_generales","brillo._calibrada.pdf"), format='pdf', bbox_inches='tight')
    plt.close(fig6)
    
    ######################################################################
    # Ajustes con imfit (exponencial)
    ######################################################################
    # creamos una carpeta donde guardar los outputs:
    if os.path.exists("exponencial"):
        shutil.rmtree("exponencial")  
    os.makedirs("exponencial") 

    # Ajuste con una exponencial solamente:
    residuo = os.path.join("exponencial", "exp_residuo.fits")
    modelo = os.path.join("exponencial", "exp_modelo.fits")
    params = os.path.join("exponencial", "exp_params.dat")
    comando = (
            f'cd {os.getcwd()} && '
            f'imfit {galaxia} -c {config_exp} --sky={sky} --gain={gain} --readnoise={ron} --save-model={modelo} --psf=mofeta.fits --mask={mask} --save-residual={residuo}  --save-params={params}'
            )
    result=subprocess.run(comando, shell=True, capture_output=True)
    
    # Leemos el fichero de parámetros para poder plotear el modelo exponencial
    # encima de la curva de brillo superficial:
    with open(params, 'r') as params:
        lines = params.readlines()
        print(lines)
        I_0_exp = [float(l.split()[1].strip()) for l in lines if "I_0" in l][0]
        h_exp = [float(l.split()[1].strip()) for l in lines if l[0]=="h"][0]
    
    # plot del modelo exponencial:
    hdu_list_model_exp = fits.open(modelo)
    model_exp = hdu_list_model_exp[0].data
    hdu_list_model_exp.close()
    model_exp_cal = -2.5*np.log10(model_exp/pixel_size**2)-zcal

    fig7, ax7 = plt.subplots(figsize=(12, 10))
    cax = ax7.imshow(model_exp_cal, cmap='inferno', origin='lower', aspect="equal", extent=[x_extent.min(), x_extent.max(), y_extent.min(), y_extent.max()], 
                     vmin=np.min(calibrada_abs), vmax=23)
    cbar = plt.colorbar(cax, ax=ax7, pad=0.13)
    cbar.set_label(r"$\mu$ [$\unit{\mag\per\arcsec\squared}$]", fontsize=30)
    cbar.ax.tick_params(labelsize=24)
    ax7.set_xlabel(r"$X$ [$\unit{\arcsec}$] ", fontsize=30)
    ax7.set_ylabel(r"$Y$ [$\unit{\arcsec}$]", fontsize=30)
    ax7.tick_params(axis='both', which='major', labelsize=24)
    
    kpcax = ax7.secondary_xaxis('top', functions=(arcsec_kpc, kpc_arcsec))
    kpcax.set_xlabel(r"$X$ [$\unit{\kilo\parsec}$] ", fontsize=30, labelpad=14)
    kpcax.tick_params(axis='both', which='major', labelsize=24)
    
    kpcay = ax7.secondary_yaxis('right', functions=(arcsec_kpc, kpc_arcsec))
    kpcay.set_ylabel(r"$Y$ [$\unit{\kilo\parsec}$]", fontsize=30, labelpad=14)
    kpcay.tick_params(axis='both', which='major', labelsize=24)
    
    fig7.savefig(os.path.join("exponencial","modelo_exp.pdf"), format='pdf', bbox_inches='tight')
    plt.close(fig7)
    
    #plot del residuo:
    #hdu_list_resid_exp = fits.open(residuo)
    #resid_exp = hdu_list_resid_exp[0].data
    #hdu_list_resid_exp.close()
    #resid_exp_cal = resid_exp
    #resid_exp_cal[resid_exp_cal>0] = -2.5*np.log10(resid_exp_cal[resid_exp_cal>0]/pixel_size**2)-zcal
    #resid_exp_cal[resid_exp_cal<0] = +2.5*np.log10(-resid_exp_cal[resid_exp_cal<0]/pixel_size**2)+zcal
    resid_exp_cal  = calibrada_abs - model_exp_cal

    fig8, ax8 = plt.subplots(figsize=(12, 10))
    cax = ax8.imshow(resid_exp_cal, cmap='PRGn', origin='lower', aspect="equal", extent=[x_extent.min(), x_extent.max(), y_extent.min(), y_extent.max()],
                     vmin=-0.5, vmax=0.5)

                     
    cbar = plt.colorbar(cax, ax=ax8, pad=0.13)
    cbar.set_label(r"$\mu$ [$\unit{\mag\per\arcsec\squared}$]", fontsize=30)
    cbar.ax.tick_params(labelsize=24)
    ax8.set_xlabel(r"$X$ [$\unit{\arcsec}$] ", fontsize=30)
    ax8.set_ylabel(r"$Y$ [$\unit{\arcsec}$]", fontsize=30)
    ax8.tick_params(axis='both', which='major', labelsize=24)
    
    kpcax = ax8.secondary_xaxis('top', functions=(arcsec_kpc, kpc_arcsec))
    kpcax.set_xlabel(r"$X$ [$\unit{\kilo\parsec}$] ", fontsize=30, labelpad=14)
    kpcax.tick_params(axis='both', which='major', labelsize=24)
    
    kpcay = ax8.secondary_yaxis('right', functions=(arcsec_kpc, kpc_arcsec))
    kpcay.set_ylabel(r"$Y$ [$\unit{\kilo\parsec}$]", fontsize=30, labelpad=14)
    kpcay.tick_params(axis='both', which='major', labelsize=24)
    
    fig8.savefig(os.path.join("exponencial","residuo_exp.pdf"), format='pdf', bbox_inches='tight')
    #plt.close(fig8)
    
    #plot del modelo exponencial con la curva de brillo superficial:
    rr = np.linspace(0,50,200)
    mu_exp = -2.5*np.log10(I_0_exp*np.exp(-rr/(h_exp*pixel_size)) /pixel_size**2)-zcal
    fig9, ax9 = plt.subplots(figsize=(14, 10))
    ax9.plot(smas_cal, mu, "o", color="purple", markersize=9, alpha=0.7, label=r"Medidas")
    ax9.plot(rr, mu_exp, "--", color="darkgreen", linewidth=3, alpha=0.7,label=r"Modelo exponencial")
    ax9.set_xlabel(r"semieje mayor [$\unit{\arcsec}$] ", fontsize=30)
    ax9.set_ylabel(r"$\mu$ [$\unit{\mag\per\arcsec\squared}$]", fontsize=30)
    ax9.tick_params(axis='both', which='major', labelsize=24)
    ax9.invert_yaxis()
    ax9.legend(fontsize=24)
    
    kpcax = ax9.secondary_xaxis('top', functions=(arcsec_kpc, kpc_arcsec))
    kpcax.set_xlabel(r"semieje mayor [$\unit{\kilo\parsec}$] ", fontsize=30, labelpad=14)
    kpcax.tick_params(axis='both', which='major', labelsize=24)
    
    fig9.savefig(os.path.join("exponencial","curva_brillo.pdf"), format='pdf', bbox_inches='tight')
    plt.close(fig9)
    ######################################################################
    # Ajustes con imfit (exponencial+sersic)
    ######################################################################
    if os.path.exists("exponencial+sersic"):
        shutil.rmtree("exponencial+sersic")  
    os.makedirs("exponencial+sersic") 
    
    # Ajuste con una exponencial + sersic:
    residuo = os.path.join("exponencial+sersic", "expser_residuo.fits")
    modelo = os.path.join("exponencial+sersic", "expser_modelo.fits")
    params = os.path.join("exponencial+sersic", "expser_params.dat")
    comando = (
            f'cd {os.getcwd()} && '
            f'imfit {galaxia} -c {config_model} --sky={sky} --gain={gain} --readnoise={ron} --save-model={modelo} --psf=mofeta.fits --mask={mask} --save-residual={residuo}  --save-params={params}'
            )
    result=subprocess.run(comando, shell=True, capture_output=True)
    
    # Leemos el fichero de parámetros para poder plotear el modelo completo
    # encima de la curva de brillo superficial:
    with open(params, 'r') as params:
        lines = params.readlines()
        print(lines)
        I_0_exp = [float(l.split()[1].strip()) for l in lines if "I_0" in l][0]
        I_0_exp_err = [float(l.split()[4].strip()) for l in lines if "I_0" in l][0]
        
        h_exp = [float(l.split()[1].strip()) for l in lines if l[0]=="h"][0]
        h_exp_err = [float(l.split()[4].strip()) for l in lines if l[0]=="h"][0]
        
        n_ser = [float(l.split()[1].strip()) for l in lines if l[0]=="n"][0]
        n_ser_err = [float(l.split()[4].strip()) for l in lines if l[0]=="n"][0]
        
        I_e_ser = [float(l.split()[1].strip()) for l in lines if "I_e" in l][0]
        I_e_ser_err = [float(l.split()[4].strip()) for l in lines if "I_e" in l][0]
        
        r_e_ser = [float(l.split()[1].strip()) for l in lines if "r_e" in l][0]
        r_e_ser_err = [float(l.split()[4].strip()) for l in lines if "r_e" in l][0]
        
    # plot del modelo exponencial + sersic:
    hdu_list_model = fits.open(modelo)
    model_expser = hdu_list_model[0].data
    hdu_list_model.close()
    model_expser_cal = -2.5*np.log10(model_expser/pixel_size**2)-zcal

    fig10, ax10 = plt.subplots(figsize=(12, 10))
    cax = ax10.imshow(model_expser_cal, cmap='inferno', origin='lower', aspect="equal", extent=[x_extent.min(), x_extent.max(), y_extent.min(), y_extent.max()], 
                     vmin=np.min(calibrada_abs), vmax=23)
    cbar = plt.colorbar(cax, ax=ax10, pad=0.13)
    cbar.set_label(r"$\mu$ [$\unit{\mag\per\arcsec\squared}$]", fontsize=30)
    cbar.ax.tick_params(labelsize=24)
    ax10.set_xlabel(r"$X$ [$\unit{\arcsec}$] ", fontsize=30)
    ax10.set_ylabel(r"$Y$ [$\unit{\arcsec}$]", fontsize=30)
    ax10.tick_params(axis='both', which='major', labelsize=24)
    
    kpcax = ax10.secondary_xaxis('top', functions=(arcsec_kpc, kpc_arcsec))
    kpcax.set_xlabel(r"$X$ [$\unit{\kilo\parsec}$] ", fontsize=30, labelpad=14)
    kpcax.tick_params(axis='both', which='major', labelsize=24)
    
    kpcay = ax10.secondary_yaxis('right', functions=(arcsec_kpc, kpc_arcsec))
    kpcay.set_ylabel(r"$Y$ [$\unit{\kilo\parsec}$]", fontsize=30, labelpad=14)
    kpcay.tick_params(axis='both', which='major', labelsize=24)
    
    fig10.savefig(os.path.join("exponencial+sersic","modelo_expser.pdf"), format='pdf', bbox_inches='tight')
    plt.close(fig10)
    
    #plot del residuo:
    # hdu_list_resid_expser = fits.open(residuo)
    # resid_expser = hdu_list_resid_expser[0].data
    # hdu_list_resid_expser.close()
    # resid_expser_cal = resid_expser
    # resid_expser_cal[resid_expser_cal>0] = -2.5*np.log10(resid_expser_cal[resid_expser_cal>0]/pixel_size**2)-zcal
    # resid_expser_cal[resid_expser_cal<0] = +2.5*np.log10(-resid_expser_cal[resid_expser_cal<0]/pixel_size**2)+zcal
    
    resid_expser_cal  = calibrada_abs - model_expser_cal
    fig11, ax11 = plt.subplots(figsize=(12, 10))
    cax = ax11.imshow(resid_expser_cal, cmap='PRGn', origin='lower', aspect="equal", extent=[x_extent.min(), x_extent.max(), y_extent.min(), y_extent.max()],
                      vmin=-0.5, vmax=0.5)

                     
    cbar = plt.colorbar(cax, ax=ax11, pad=0.13)
    cbar.set_label(r"$\mu$ [$\unit{\mag\per\arcsec\squared}$]", fontsize=30)
    cbar.ax.tick_params(labelsize=24)
    ax11.set_xlabel(r"$X$ [$\unit{\arcsec}$] ", fontsize=30)
    ax11.set_ylabel(r"$Y$ [$\unit{\arcsec}$]", fontsize=30)
    ax11.tick_params(axis='both', which='major', labelsize=24)
    
    kpcax = ax11.secondary_xaxis('top', functions=(arcsec_kpc, kpc_arcsec))
    kpcax.set_xlabel(r"$X$ [$\unit{\kilo\parsec}$] ", fontsize=30, labelpad=14)
    kpcax.tick_params(axis='both', which='major', labelsize=24)
    
    kpcay = ax11.secondary_yaxis('right', functions=(arcsec_kpc, kpc_arcsec))
    kpcay.set_ylabel(r"$Y$ [$\unit{\kilo\parsec}$]", fontsize=30, labelpad=14)
    kpcay.tick_params(axis='both', which='major', labelsize=24)
    
    fig11.savefig(os.path.join("exponencial+sersic","residuo_expser.pdf"), format='pdf', bbox_inches='tight')
    #plt.close(fig11)
    
    #plot del modelo exponencial + sersic con la curva de brillo superficial:
    rr = np.linspace(0,50,200)
    mu_exp = -2.5*np.log10(I_0_exp*np.exp(-rr/(h_exp*pixel_size)) /pixel_size**2)-zcal
    
    b_n = 1.9992*n_ser - 0.3271
    b_n_err = 1.9992*n_ser_err
    mu_sersic = -2.5*np.log10(I_e_ser*np.exp(-b_n*((rr/(r_e_ser*pixel_size))**(1/n_ser)-1)) /pixel_size**2)-zcal
    
    mu_total = -2.5*np.log10(((I_e_ser*np.exp(-b_n*((rr/(r_e_ser*pixel_size))**(1/n_ser)-1))) + I_0_exp*np.exp(-rr/(h_exp*pixel_size))) /pixel_size**2)-zcal 
    
    fig12, ax12 = plt.subplots(figsize=(15, 10))
    ax12.plot(smas_cal, mu, "o", color="purple", markersize=9, alpha=0.7, label=r"Medidas")
    ax12.plot(rr, mu_exp, linestyle="dashdot", color="darkgreen", linewidth=3, alpha=0.8,label=r"Perfil exponencial")
    ax12.plot(rr, mu_sersic, linestyle="dotted", color="darkorange", linewidth=3, alpha=0.8,label=r"Perfil de Sérsic")
    ax12.plot(rr, mu_total, "--", color="black", linewidth=3, alpha=0.9,label=r"Total")
    ax12.set_xlabel(r"semieje mayor [$\unit{\arcsec}$] ", fontsize=30)
    ax12.set_ylabel(r"$\mu$ [$\unit{\mag\per\arcsec\squared}$]", fontsize=30)
    ax12.tick_params(axis='both', which='major', labelsize=24)
    ax12.invert_yaxis()
    ax12.legend(fontsize=24, loc=3)

    # Inset plot
    inset = inset_axes(ax12, width="50%", height="40%", loc="upper right")  # Adjust size and location
    inset.plot(smas_cal, mu, "o", color="purple", markersize=9, alpha=0.7, label=r"Medidas")
    inset.plot(rr, mu_exp, linestyle="dashdot", color="darkgreen", linewidth=3, alpha=0.8,label=r"Perfil exponencial")
    inset.plot(rr, mu_sersic, linestyle="dotted", color="darkorange", linewidth=3, alpha=0.8,label=r"Perfil de Sérsic")
    inset.plot(rr, mu_total, "--", color="black", linewidth=3, alpha=0.9,label=r"Total")
    inset.set_xlim(0, 10)  # Set the zoomed-in range
    inset.set_ylim(16, 13.5)
    inset.set_xlabel(r"semieje mayor [$\unit{\arcsec}$] ", fontsize=20)
    inset.set_ylabel(r"$\mu$ [$\unit{\mag\per\arcsec\squared}$]", fontsize=20)
    inset.tick_params(axis='both', which='major', labelsize=18)
    
    kpcax = ax12.secondary_xaxis('top', functions=(arcsec_kpc, kpc_arcsec))
    kpcax.set_xlabel(r"semieje mayor [$\unit{\kilo\parsec}$] ", fontsize=30, labelpad=14)
    kpcax.tick_params(axis='both', which='major', labelsize=24)
    
    fig12.savefig(os.path.join("exponencial+sersic","curva_brillo.pdf"), format='pdf', bbox_inches='tight')
    plt.close(fig12)
    
    ######################################################################
    # Cocientes B/T, D/T
    ######################################################################
    # Luminosidad total disco:
    L_disc = 2*np.pi*I_0_exp*h_exp**2
    L_disc_min = 2*np.pi*(I_0_exp-I_0_exp_err)*(h_exp-h_exp_err)**2
    L_disc_max = 2*np.pi*(I_0_exp+I_0_exp_err)*(h_exp+h_exp_err)**2
    L_disc_err = abs(L_disc_max - L_disc_min)/2
    
    L_bulbo = 2*np.pi*I_e_ser*r_e_ser**2*n_ser*np.exp(b_n)*gamma(2*n_ser)/(b_n**(2*n_ser))
    L_bulbo_min = 2*np.pi*(I_e_ser-I_e_ser_err)*(r_e_ser-r_e_ser_err)**2*(n_ser-n_ser_err)*np.exp(b_n-b_n_err)*gamma(2*(n_ser-n_ser_err))/((b_n-b_n_err)**(2*(n_ser-n_ser_err)))
    L_bulbo_max = 2*np.pi*(I_e_ser+I_e_ser_err)*(r_e_ser+r_e_ser_err)**2*(n_ser+n_ser_err)*np.exp(b_n+b_n_err)*gamma(2*(n_ser+n_ser_err))/((b_n+b_n_err)**(2*(n_ser+n_ser_err)))
    L_bulbo_err = abs(L_bulbo_min - L_bulbo_max)/2
    
    L_tot = L_disc + L_bulbo
    L_tot_err = L_disc_err + L_bulbo_err
    
    D_T = L_disc / (L_disc+L_bulbo)
    D_T_err = D_T*np.sqrt((L_disc_err/L_disc)**2 + (L_tot_err/L_tot)**2)
    
    B_T = L_bulbo / (L_disc+L_bulbo)
    B_T_err = B_T*np.sqrt((L_bulbo_err/L_bulbo)**2 + (L_tot_err/L_tot)**2)
    
    B_D = L_bulbo / L_disc
    B_D_err = B_D*np.sqrt((L_disc_err/L_disc)**2 + (L_bulbo_err/L_bulbo)**2)
    
    print(f"D_T = {D_T} +- {D_T_err}")
    print(f"B_T = {B_T} +- {B_T_err}")
    print(f"B_D = {B_D} +- {B_D_err}")
    
