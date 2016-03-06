import wx
import os
import os.path
import platform
import sys
from wx.lib.embeddedimage import PyEmbeddedImage
sys.path.append("/home/balto/VirtualBox VMs/shared/InteractiveROSETTA/dist/InteractiveROSETTA/eggs")
import poster
#from poster.encode import multipart_encode
#from poster.streaminghttp import register_openers
import urllib2
import glob
import socket

poster.streaminghttp.register_openers()
if (platform.system() == "Windows"):
    # The reason for the slowness at Rosetta launch is because it's loading the Dunbrack rotamers
    # from text files.  It attempts to write out a binary file of these rotamers after loading them
    # from text once, but for some reason on Windows it doesn't write out the binary file properly
    # so the next time it thinks the binary file is "outdated" and writes a new one every time!
    # Use the server for Windows to get around this inconvenience
    # We also have to go through the User's AppData directory and delete all these temporary
    # binary files because each one uses 70MB, so disk space will be lost quickly if they aren't
    # cleaned up
    useServer = True
    homedir = os.path.expanduser("~")
    binaries = glob.glob(homedir + "\\AppData\\Local\\Temp\\dun10_binary*")
    for binary in binaries:
	try:
	    os.remove(binary)
	except:
	    pass
else:
    useServer = False

serverName = [""]
serverTimeout = 60*10
primaryRender = ["cartoon", "ribbon"]

def getServerName():
    return serverName[0]

def setServerName(name):
    serverName[0] = name

def getPrimaryRender():
    return primaryRender[0]

def setPrimaryRender(renderType):
    if (renderType == "cartoon"):
	primaryRender[0] = "cartoon"
	primaryRender[1] = "ribbon"
    else:
	primaryRender[0] = "ribbon"
	primaryRender[1] = "cartoon"

colors = ["0xFF0000", "0x0000FF", "0x00FF00", "0xFFFF00", "0xFF00FF", "0x00FFFF", "0xFF9900", "0x6666FF",
	  "0x800000", "0x000099", "0x006600", "0xFF99FF", "0x33CCCC", "0x663300", "0x00CCFF", "0xCCCC00",
	  "0x99FF33", "0x99CCFF", "0x660066", "0xFF3399"]
def getChainColor(row):
    return colors[row % len(colors)]

# Useful scoretype translations for the user
scoretypes = {}
scoretypes["fa_atr"] = "VDW Attraction Energy"
scoretypes["fa_rep"] = "VDW Repulsion Energy"
scoretypes["fa_sol"] = "Solvation Energy"
scoretypes["fa_intra_rep"] = "Intraresidue VDW Repulsion Energy"
scoretypes["fa_intra_atr"] = "Intraresidue VDW Attraction Energy"
scoretypes["fa_elec"] = "Electrostatics Energy"
scoretypes["pro_close"] = "Proline Closure Energy"
scoretypes["hbond_sr_bb"] = "Hydrogen Bond Energy (Short Range BB)"
scoretypes["hbond_lr_bb"] = "Hydrogen Bond Energy (Long Range BB)"
scoretypes["hbond_bb_sc"] = "Hydrogen Bond Energy (BB to Sidechain)"
scoretypes["hbond_sc"] = "Hydrogen Bond Energy (Sidechain)"
scoretypes["dslf_fa13"] = "Disulfide Bond Energy"
scoretypes["rama"] = "Ramachandran Energy"
scoretypes["omega"] = "Omega Dihedral Energy"
scoretypes["fa_dun"] = "Dunbrack Rotamer Energy"
scoretypes["p_aa_pp"] = "Probability of AA Given Phi-Psi Energy"
scoretypes["ref"] = "Reference Energy"
scoretypes["total_score"] = "Total Energy"

icon = PyEmbeddedImage(
    "iVBORw0KGgoAAAANSUhEUgAAAEAAAABACAYAAACqaXHeAAAABHNCSVQICAgIfAhkiAAAAAlw"
    "SFlzAAALEgAACxIB0t1+/AAAABh0RVh0VVJMAGh0dHA6Ly93d3cucHltb2wub3JnmI83TgAA"
    "ABx0RVh0U29mdHdhcmUAQWRvYmUgRmlyZXdvcmtzIENTNui8sowAAA5bcHJWV3ic1VtNjCPF"
    "FW57drvHu3HX2Duxu2fbO82u+2fdPYu7y0sOHFAue17lFIkDyiGOhITgxDXnvXEAJEuEAwIu"
    "XKKwt4icEELRHCIOwA0QewKhKBdEpCid91512+7/8uwso5SnZ+yxXV99r169n6rXf//v3/6p"
    "vKK8kmA7hUdyenpKz07pOj1dweN0tVrRs9UqWdGf1Wq5Wp7Cr9VyueL4YrlMlsvTJT1bLvmS"
    "J0v4P4enfKmslgo8TTi9wXnC+SnnK04vOFe4suLwAQWeKvQrUZRTRVkpyhJeK/S6fBmGkbz8"
    "8ssPYeSKzAXtQV1f53n1+/3krbfeSr766qvku+++S7799tvk3XffTSzLeqx+r169mty/f//B"
    "Dnxf+Dn46rqevPPOO8k333yTfP/998nXX39N/K9du/ZY/e7v7yfPPvvsCz8T35+7/b/idy/1"
    "hvrQsibHgR3gz7FnHRoD1lM7u+FrmpbcvXv3V+ncybaO2mODQ2Ny7BH8zPZsyzhkg97l7m74"
    "ys78NXVoTWaCNfya4Q/wT18fB9ZkoGly4GXdlcHXBraVCj0Vvh2kUqC/1lCVw1dS3leuXEle"
    "ffVViQGoDKgHAjMOQh4uNo2HYSoQuCbWUAJ8Z/5d3bYAHIU9C4Mw5Hn8UEwCXJatX2rVhIcP"
    "HyaPHj1Kfvjhh+Szzz5rGUCXDVDGyDEMcsw3Q4BBcJICTo03aJHBFu+/KC361/NSDQurofPj"
    "SKciCHp6Q58vvviitP7rtkcKVxB6dQsXICIcr6v3WkQgg9/VLIvYhEEr9NYghAgsa099THyt"
    "a1lCryS4b+ZALBPL0hoUUQZfPxam5SQHbmJrGcKCvmgF9ZPQjj8cCNHnevYZA8OtMeb7LZMg"
    "VsNwcHb8oVD7HPxC4Gvt+JwsY+1KbMO3hmTpFvl5N000sN3LmtY2BTQJJAHrTPiHhrB0+T5N"
    "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"
    "ewL4elSn9EG3T17P8sqw0T4T97M/c515j8p8rrbwf99HKZ9QpjwvMdYxfjKW0Z7Xhp9SL8IY"
    "N/iPABvTvhBzbfd/H3Nyj/KY//l/IvMo9fvd/7Myn6tj/s+5HTv0fpJ1LfXxKX2Dv4jLPLZV"
    "+DG7Zxi25P0652HGcOJi57Q1e534M/coj5WDf2vxIW0nbcqe2cj/ozKf8y7IflvWKX1H3866"
    "Yo/RWEXcTK/n1/3Z+8GacMKW6pVh1IO5pPs35/LRNxjP9+dGefUw2kDfi0wbEz/znpW597VL"
    "aGm9QD2+9L9SSimllFJKKaWUUkpdTTsRERERERERERERERERERERERERERERERERERERERER"
    "ERERERERERERERERERERERERERERERERERERERERERERERERERERERERERERERERERERkTvk"
    "H4eXjmrZO46cAAABU21rQlT6zsr+AH+lhQAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAeJzt1uFpg2AUhlEHcREH"
    "cRAXcRAHcREHsbyBC7emIf+KCeeBQ5tP++tNbM5TkiRJkiRJkiRJkiRJkiRJkiRJH9FxHOe+"
    "70/nOcu1d/e/uk/3b13XcxzHc5qmx8/sGP0s99S9dRbLsjxexzAMf76HdO+yY5V9s2F2rc37"
    "PbV/1Te//o3uX7bre1Y565/lep19+8bZv7pe0/3Lc77vX//X53l+2j/X7P99Zdt67tfv27b9"
    "+sz357/9v6/6Htf3q/dArtV3+5xF1Z8d12uSJEmSJEmSJEn69wYAAAAAAAAAAAAAAAAAAAAA"
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAPhAPwr5rLhS2ipmAAACuG1rQlT6zsr+AH++HgAA"
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    "AAAAAAAAAAAAAAAAeJzt3T9OG0EUB2CQXVFT0NHRcQZXHCNU1JY4gdsUlmiDxB2QUtByhFSR"
    "KJDcULpBQkq5mYksJzEmxo5330P+iq8C27PvN/zZmdmZvaZp9lp0fNa//FR8Lm6L78W0aBZM"
    "Z1+7nX1vfU19bZttY/v2S2aD4kvxuCTndT3O3qu+577+kNZhyWZUTLaQ+Vsms8841A/SOCpZ"
    "jIuXFnNf9DL7zCP9IEy/1H5YPHeY+6LnWRv6+kGnTku9vwXmvqi25VQf6MRFqfOPBJkvqm26"
    "0Ada0yu1vUmQ8yq1jT39YKsOSj2/Jsj2vWpbD/SBrWV/nyDTdd3rA/+t98F+7pf9HvC3YHMf"
    "4e/9Kjfy38hFguy2xX3Beuq9dMZ7vE3VazE+8D51LC3T2M621GsyTrjaMEFWbRnK/5/qfErk"
    "eH7b6rWZM3rbOEFGbRvLf6k6p97lHG6Ueo3WD7w2SpBNV0by/0tdVzVJkEtX6rVaS/bbIEEm"
    "XRvIf+46QR5du5b/3FOCPLr2JP9fThJkEeVEH2jOE+QQ5Vz+zVWCHKJcyb+5S5BDlDv5Nw8J"
    "cojyIP9m2bOYu2Iq/yY6g2jR9Y8WXf9o0fWPFl3/aNH1jxZd/2jR9Y/m/7/4DCK5/4vPIJLx"
    "n/gMIhn/jc8gkvmf+Awimf+NzyCa9R+7zfqv3TZIkEfXBvKfs/6bUYJcujKS/Sue/8Lzn7vN"
    "89/Y/2G32f8F+z9h/zfs/7jb7P+K/Z+x/zvOf6By/gvOf8L5b1TOf6Ry/iuV85/50/GZ898B"
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgNx+AqUyY7s6wM9HAAARs21rQlT6zsr+"
    "AH/FhwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    "AAAAAAAAAAAAAAAAAAAAAAAAeJztnQd0lNXWhieN0JRID4QEkABC6KFEWmihiCAqVWAJSEdR"
    "il5BFBHwiiUIAibCFSwREX6keKVJuQapkiCg1ARIwdAJLSHl/d99mLCAkDoTE9baWetZaTPf"
    "N3P2Prudfc5YAFjsjIPFYnEiRYg/+RdZS86QFAIrKda/rbU+xt/6HHmuXMPer0vJOxyssvMg"
    "fUgw2UViyFWSZJV36l3yT7X+Lcn6mBjrc4Kt1/CwXlN1oeDiSNGUIJ3IGhJPku+S9/0yz4jU"
    "ux6fZL1GvPWanaz3cFQ9KDjwqyTpScLvs+t3cHZ2hqurq/meDR3IiBTrPXpa72lRPchXinL4"
    "25EdmcmtRIkS6N69O0aOHIlu3bqhbNmycHR0tEUPYL2n3Luo6kC+4M5hn0USM5NTuXLl8Omn"
    "n+KPP/7AgQMH8MsvvyAwMBA1atSAg4ODrTqQaH0N7qoD/ygtONzbs5KPzPGXXnoJJ0+exJUr"
    "V3Dx4kWcPn0av/76K1577TUUKlTIVvmnIa+lhepAniP52DPkdHbkUrx4ccyfPx+XL1+GfN26"
    "dQvnz5/Hzp07MXPmTDz66KP2kj+sr0lem5PqQZ7gzGEdZrmdn2VLJjK/33//ffz999+4du2a"
    "mf9//fUX1qxZg7lz56J585a0EXaTP6yvTV6js+qA3ZFxTcipTJ566ils3boVx48fN/7/p59+"
    "ws8//4zw8ANYvjwE/v4+Obqe+JQiRYpk9hh5jSNU/nalu+V2Dp7jOVm4cGH06dMHCxcuxDff"
    "fIP169fj1KnTuHnzPK5fP44ZM0L4GB/Ggk5ZXsvDwwO9evXCmDFj0LlzZ5QuXTozO9BddcAu"
    "+HEYo22xyxLnV61aFW3btsWPP/6IGzcuIjk5CqtXX0W7dgmoWDEIbm49mSc+jWLFGnGOF7tr"
    "vjugaNGieOKJJ0zecOjQIRw7dgx79+7F9OnTUaFChYzuK6/ZT3XAJipy+Hba0T8z/38KW7b8"
    "QnsQj4YNwTl9BtWr70bNmn+Y748/vhJlyoyAk1M96kNL9O07AFOnTsbixYvx559/4urVq0hI"
    "SKD9uGl0YeDAgXBxccnofvLaK6oO5JpAe8o+DV/fZ2kP9sPd/RRlvgN16pzg35LQoMEV/nyK"
    "/9sEL68Q9O9/HEeOnGHecB3nzp3lz0eN7OVL4smIiAgTX5YqVSqz+wWq/HOFvyWL2k7ucaTP"
    "n0ib/j80aXKcMeBVcgstWtxC48bnULnybjRtuhXLlgFJSUbc1IFEyv8Y4uLiEB8fj5iYGBw8"
    "eBDvvPMOfYdbZveT9yDvxaJ6kG2KWbKo6dpK6dIDOdd/R4cOl9CjB9C1K9C+PfDkkzdoE8Lh"
    "47Oc/v08tm8HEhMvITY2EqGhvxm/f/ToUeMLVq5ciU6dOtFXZBk7ylpicZV/tpA11r55KXsH"
    "B0eULDkY9esfo8zP4tlnkxjPJ6FlS5n/Fyn7A4wLvqMebKR/D8WSJQs5z8dj2LDh+PDDD7F0"
    "6VJ8/vnnfN6zjBeLZfe+8p50/ThrZD3lQN7K30K5+aFatXVG1o0aRcHP7xK/x/P3SPqFHZR/"
    "MB/3HB/XBJUqVWIsWMzk/V5eXnxcI8aJj2cW9z0IeU+6VpQ5Mj8C8lL2aTg5FWecP4H523/p"
    "C7bB03MP5X6IecDvlO0akwNYLHZbH0gjQG1ApkjddPU/IX/B3X08nnsuDIMH70KtWitMzO/p"
    "GWRk7+Linhf3lPemteGMKW+5q8YrcVX58uUpJ3dTx7OvLB5D3bobmdffxLFjNzBmTBjt/Fy4"
    "ujbkfR/JK52T91Ze5f9AxC72slj7d8qUKYNXXnnF1F2kdvv666+jdu3adpOFi0tLdOsWiU2b"
    "bmHnzmQMGXIG5cq9ndc2R95bb4v6gIzk/wVJkjjrvffeMzlWdHS0WbPbsWMHPv74Y8Zs1ewi"
    "i2bNumD+/K04fPg6vvwyjrl/GIoW7WyXa0vPSZMmTVCvXj088sg9tkT6CReq/B+I9NfucXR0"
    "TPb19cWePXtw6dIlU2s5d+6c0YEVK1aYPi7bZeSA0aNfwsGD23HzZjxCQtaiY8d/wdm5jM3X"
    "rlOnDqZMmYJFixZh3rx5GDVqlPFf1v9LL+key+1e4vwe74JGKw7JGco/pUuXLqa2Lv0aSUlJ"
    "pndHfl+2bBnz8QH0z7nv4ZQ4onr16nj11VcRFraXuhVHOQWiTZtGZq3HFtlLn+G4ceOwevVq"
    "s84ohISEmDqBtUaUts+glco/HZM4JNdIqvTmbdiwwdRYpXdD6m3r168z9TaJB8QGuLjkXFZV"
    "qlTB008/jf79+zPufw5Dhw41fmbEiBHM+SrbPPelJvDRRx9h8+bNpsdo+/btWLJkCYYPH258"
    "gaenZ6qzs7O8x0kq/3Ss5JDcSpujIhvp1ZQxlH6NNWtWIzLyJHUiij8vQ7t27XIsG1mrGz16"
    "tJn7L7/8svm9R48epkdE1nelLmiL/CVXEdsvPSabNm1iXPElJk2aZPyA6O5XX30lunaLdmKl"
    "yj8dcZa7eveld69r165488038cUXX5i1t0uXztMfxNI2ROLJJ1fxcR3JI7TbRSi79LW4tP5e"
    "6QMTGU+YMAFvvPGGWbMR5GfJMQYNGoSAgAA89thjNslfeoM6duyIOXPmYMGCBeb6ogP79u3D"
    "iRMnTA/yli1bUnr16hXnKI1n+T/mBQnZd3PP/hzxmUWKFDZzffPmjUhIOE2bcBFt26Yyxj7E"
    "fH0hvLwWkPkoVWoo47e0XowSfO4T8PDwgZ9fIyNbsfmTJ0829ftPPvkEs2bNMv0bol/Dhg0z"
    "dkByC1t7wmV/ibzed99919xD7NfZs2dNHCu9p/v370+dOXNmqpeXV36Pd0Ejk5jNCT17tsPs"
    "2TvRunUq/WgcY7hQ2uww1Kr1F2rU2E3fvow6MQkVK0429Vtf372YOHEJfvghCMHB8zjHBxv5"
    "BwcHG3sic1TkI38Tn/D888/zWrXssSfArBOIvondDwsLw/Xr15GYmIgLFy4YGyD7Eehv8nu8"
    "CxqZjqmTU2G4uQ2Ct/c+1K79Jxo0OIrGjS+jYcPL/DmGMv+Nsv+KseEO2vabnOfXaGtjkZx8"
    "BWfORBnbL35fxl56f8UOTJ061fxd9gVJPGiP+Z9G/fr1TSwgMYzY/lOnTpk4Vn6fOHGi9Izk"
    "93gXNNLZ//Q64E6bH4RmzWLpAy7TziaiVatktGiRQL04iNKll2Ls2F3Yti0FK1aAedgFHD58"
    "kL5jHef4GAwePBhvv/22kXua/xfZ9+vXDx06dKD8vW3OAdOQvQVyP8kHpOdQ4sF169aJ7qW2"
    "bNkylb4tv8e7oHFP/Pdg3+pBO/89ZX6Z/joZzzwDyg1o2TIJdeueQoUKP8LffyXn1xnG9tf4"
    "mL0YP34Kc4k+jMs6GRsvtl5if0FkL73B0ssruYDY7GrVPOxiA7y9vU0MKPmq6Jr0icmeo969"
    "e6e0atUqjr4mv8e7oHEn/8t4/pdmTPcZmje/iC5dbpBbjAeS6AeuUP4RZj3f0/NTVK48j7HA"
    "LPqLfihb1ou5hDNcXJxNDii1A5nvIu+ePXuaeoDkmtOmSU4wljo1FEWL1su13CXXkPpv3759"
    "sXbtWhP7b9y40cQcshdtwIABt6iHK5kn5Pd4FzTu1H8yGltHR1fmaP3h47OTPiDK+AGJAerW"
    "jUXNmmGcc6tQsuQQzt+0PO7eeSzzumTJksbPy1qSIP38Ivt//3sCbcUc5pXHqUM/MOb0Tvf8"
    "rJBr+/v7Gzsj8YTsPdu9ezdCQ0PxwQcfiH1JHT9+/DXqwCTmtvk93gUNU//Nygc4OFTknF5E"
    "GYaRbZT7Tsb/O/nzBtr/GXB19c62vCRfF12QtYDRo7+mPokdAWO3C9SBObQDjagHNTmnxSdk"
    "3gsivUC06yaXFB8jNUXxL5JfTJs2zew79/X1TeH/zgwaNKgVf8/v8S5omPUfy+01kgxkL31b"
    "PpRNOPOnBMo+nLn/cri7z2Y8/SJlUClXNrtKFT/TB1avHmhPwJwinnqxkX5kKapWXc77zTd7"
    "QxwdZb+orOeVJZWY55WlrS9Jn+Ru1n3ElnB+m5qfyFxiTYkzJA6UNYAWLVok0/fsIUXU/qfj"
    "zvpvRnJydXVk/NQXISFxzN1BHwvO13D69xf4f1vO9XCmfEfyWpGmB7BGjb2oXn0XZRpNvYig"
    "nu2ljnzL/OJl6tkk2p9vKfNNCAhYwPn9CgID36PMp+CFF/qb3G727Nn47LPPzHex+6IPEmME"
    "BAQkMf5YSNvvQB+R3+Nd0Lin/+NBSK/t/PmBfOhNyFd0dBKmT/8Z5cs3sDlel/px2bLjKftQ"
    "yjycehCNpk1v0GbH0y5E0Q5spR6sYOwYjTFjgFGjgODgOEREHEJi4hXs33/Q7AscO3asWfuV"
    "dR+p/0jPgsT/4hMYa6Yw9+vt5+fnQPJ7vAsi9/R/3U+jRr5mLSgh4bb8ExL+5jybRhuc6f6L"
    "bMrfgTagG+W/mflkDPPFm8wtU9Cmze09AZ6evzEuWEWZX+VrAP7zH2D9euD69ZtITk5EePh+"
    "6sRI4/OlxiSxn8heag2yJixrDIwNE5gblNfcP0Me2P8pcZr008j6WWjoduZU4aYvaNu2rfQH"
    "srZum+xvy9+ZuUU/zv2DaN/+IuP3FMZsQLt2qfTbiSbelPk/btwpzJgB5nKJGD48grJehLlz"
    "p5r6UkBAR1NDEL8v/l/qf2IPRPaSZ3p4eGj/Z+Y8sP9b1lWDgoIQGxtrzm+RfVeyrrJ48RKO"
    "eSc71Wwd6V+aMq7ciubNI9C581XagCQj+4YNz6NWrQO0Ad8yJgwkcxkHvA83twGoVOlxxgXF"
    "TO7XmMGj1BEl92OOhyFDhpjYT2JAiQu6dOkS4Obmpr1fmZNu/4f0Tsi+K+kFSk1NxY0bN3D4"
    "8GFzjofUVHKwBycLG1CccpyK2rW30O9LDBCJBg0u0CZIfWEPypR5jY+THrFCuH9vgJwt9+KL"
    "Lxo/3759e7Ru3drMebEF1jrQgVWrVhW9fPlyfo9vQSfd/i/J0Xft2gX5Sk5ONuf5REZGmnOc"
    "JM+2X2+4C33JdOrAIcYU25lbrmcesB3e3htog95EoUIVM3yurPtK/vf1118bWyXrytL3Q5mb"
    "dd+YmJi+1FuHAjC+DwP37P+U+f3WW2+ZfjDpBYyKijLrqFJbb9q0qZ1kL+d7VGFOuZXxXyrz"
    "vXi4u4eZ2kLx4q2suX/mz5eeddFHif8l95OcT34+efLkroiIiOIFYFwfJvwtd+3/lpq6rNXK"
    "WprYfTmLQ/os7HiGG+O0gcznjuP3369SbkCfPqd43145uob0Lcnar+wRlNfMXD+RNqDNiRMn"
    "LLRd+T2mDxv3nP+Qdg6LjK2cxWOHMzzvIOu1wcFBzOeiedtL5ApWrPi/u3u3c3vdQMb/+T2O"
    "Dyt2P//lQUh+IbUbsSv79u03Z31ERR3CO++MYlxn07X1/Bfbsfn8p4yQvFHidukDkf1FsudE"
    "8sodO3bhu++WmfMBbbi+nv9kP3J9/ltW8q9Zs6bp15YzYqU/T853kfrC999/b3qGc1lb0PPf"
    "7E+uzn/MCh8fH2zbts3sMZPagvRpSp4hZ0ZK3S4X19TzH/OGHJ//mh0k7pP9xVJXlLqC9Gof"
    "OXLE5O3Sy5GLea/nv+YdOTr/Obs0a9bM7NNLiwEkx5S6bQ7rSnr+8z9Hts5/zy6SR0q/pvQF"
    "Sj25efPmOa0p6/nv/zzZ+vyHnCC5gOzdyEHMp5//kL9k6/Nf8gj9/JcCgiUbn/9kJ/Tznwou"
    "d3/+m3yOo70//22tRT//7WHg/s9/DCK7Sawl+5//GGt9TpBFP//xYUY//1VRFEVRFEVRFEVR"
    "FEVRFEVRFEVRFEVRFEVRFEVRFEVRFEVRFEVRFEVRFEVRFEVRFEVRFEVRFEVRFEVRFEVRFEVR"
    "FEVRFEVRFEVRFEVRFEVRFEVRFEVRFEVRFEVRFEVRFEVRFEVRlILM/wMYH3DTLP+BLQAAKhdt"
    "a0JU+s7K/gB/1PAAAAABAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAHic7X0ruOwo1vaSSCwSicQikUgkFhmJxCIjkVgk"
    "EhmJjYyMjI0smX9R+5zunp7p+dT/1Ihac+k+VXvXCbAu77suVObnfTaeANqzkS3G10Zgh6PD"
    "AnBdxQVrAN+FfsPzYh3ggQoQAbYKG9CeJMF33ZPZsYTB8c18c/zxQ28AlZvdQSvVcTO2vmxP"
    "FRTgeJ1A4SjpMPBhua8rP/cJEqDcVCykX40DrzeBuHNcndvez5heQmwxKfxDEfOV0g8PK9Rr"
    "2yjuRnlOIjj1lmRQQ8xfORbI0j5PBjAmbKs0uI9JbSv+7utukHfu20cXj3LFsPiNmeABPFGq"
    "g3EJD9EUCSuvl7KFSJN9DPqhrsFlobcdf3GPua5+foJbKS6jNWODiTYs1vq4xcDBgm0Onh0E"
    "dU+g+O+oOXBc+NP9PC8bDy8/vPy3uE7EOhKek03CmwVwKbYVIBX2xJwtHNUeMnDAJw+HdUtx"
    "YAK+tM1ft+Da5sAf1S+4mfs2/DQdPH4AhQu0Hjc3U+obgcfhTt3VQlHX4dbt8+unqJR1TeD3"
    "e4+O+zXIJS5Cpk7JigsYazoYCWubTsC8bYE52A/85wIqp3WBVcV8MqiG2SU70e8RgZurHbhd"
    "RuFh15IpzwuqUkUlSFdjME1nA8Y+u/gpL3RpaJNmmPXVCdG4WIY+ysocqBLLRcvF8uMpFZbU"
    "PA8s6Tb2czTF4cB/1jWbeuBi8D+kokof8OD2XBs8GU8cTSVPIyg35DbgOqcWPQmdqur904sH"
    "WUGj98KDSA22qwiQTKBzNpvOA02DWOrI+UJjWJ0mx5hKvRN0BGW7Lsr2EvyozwkzLhhqZSiU"
    "zz/UPD+dLTHpJHCdTwE9AP1/eBQaEowL/9r9CR9dPEp0wqG3VmebmmB8SSw85LiVfeBG8w5R"
    "al3QbyVbUGHR/QGINv0YWBJZv8084ReqPxCoWW9oAIBGnhf8MDY34YGtHzZKRvGXR1vwhQV3"
    "dimazzc/LBzkQHeOCo0Gbk3gx6bdE23MBcprPj/16MlM2mrvD7MVPYDdD9old4NaiGl6RlR4"
    "BoEQ9IQkEYGva1D2OJtFt5Bt8vgJakFPmfHU1/regKueHD5+/pKG5dzg2IaRugbpQjn6teIJ"
    "hgvWpAI4Va2rSxwOQ8N2tGpi6w9MC+jl50O8Au+Aea8FoQvnHo07pG0XagtQLtQFIJf44+9E"
    "a/EVwup3/qFV/0XCwoAz9NyowZSRlZI4eOtVwIVKyvy5cxKPoxKJnlyEswgO6Mmfjis7Bn0H"
    "BHOtGEYQ4x1RKB5LSa3u96ZY3ZuExqgKuTELy/r+K0uP+qjoZFiMH107SsSjju9jCIh4JJ2n"
    "RNHXt94PEJ6iE1hgadceIOyo69EQQGzMj/tybrBtJIGoxl7XOc6E73pCR8+eoFE9FcZuZhDk"
    "a4RE6vasZTsKPKj9+BZh0/w+LLXiop6basbva4cwQp9bcCj14iS/HQC6h8egkdv2zHD9NAxu"
    "yxnLcWCUWMaT+Qn6ds+19ugY2S549UhujPuNb3KfSr6AzzWs8cHg/0jgHHWpifHq64eXjwtm"
    "4KcWDO3X12HsGJWGiVtaFxk6PjzHTUBKoznzAv0CrOIk03FdFQGhAH09SIUWDGsE0P4zxsoY"
    "uuOv+emyunS/UZM9f4IBLAk3xscGtd+7/ezq53MNxD6Q46Iz+Lbv3tw2W6bRZ5WolwxSTI3Y"
    "jaqo+RGtPxe3KAyNJnfdLjdDI35CewiCXa/TCtfil1XUVwKyDDeZ0jF/amt+gmWUY0e7v3IW"
    "y8f5H9DjRNguGxI99MtLtNzu6wjFQN1X3cexTRID+zDlgJAD4/vt6OS8MM5cBtryeH+Q8652"
    "z3HfTlqiCz4jBMYNg4SM4EJFlwmZpSmVgromedhBfXTlP0L76gtZ7G0owldJcOGBybHygPEL"
    "uHy9Mpcr6P3gXDK39iDt3imQbNw4t9Z0bBgFHMFAWi5CvYCj7xgElWXxhYuNg1JT3/SBxoNt"
    "PmSYSYHp/mz+9PInTg1hhmTEokczuSWNhrwjqyk/6LzPJAUBcx8c3wkDXzU9E7LtWRzHQlIj"
    "LWsicUdQLdBlEv4i52atwQjC4SXWqS3PkzMeN+rQ5MzIONRNOZkZgc+KGYosG6zo5F8qbjtI"
    "gsH6xkUWQsaxhh3WY2y/fvjO7rHnDcudW4OOL3Nhn2e4SRUXRQgy5Sx6A9Ix2hd0gRs6kmtM"
    "xtPnzsEGoc3tHMiZCA/lo4tHKeYc1HsSN8pv8MvFbmSo+KTot/DhlXtAcvVQmD4QxmvCd4xr"
    "172+oQsjuA9rWBdmeZES1kXH95rIQanNQsI5wnVNELDb3jRQPblfBNNskpDGZ1ePrtiH3U6V"
    "FNUjll9umYdH76RwA3ALLFqFHhL/VXWbNsiT98NWppvTsLjlMEVLkTcqfLf9GF2ve538NzVG"
    "XOnUtrv6elHYFaB6IeGCxwcJdRVIgD7u//OmdXCastr29VTZo7tvM1ApiPi0W+Be1Tbj1trz"
    "42AgLZpkJhLhKj22JcTAymZZkjy/XpKD2LdgXzadqN/IfGgduMzrBTPYoT6AhDIgGVC6EPpx"
    "/9c3BxXPjrML/dUO/CxOc75qu0aZPUK1ivxgC6jtgbOVQ6fy9gRpjlWSKQFS6ZCPQEzF3wbS"
    "roSL/4kdArfHp21iPDITRkiTUnGwshzDuUa9HuXj+PdYHLppjeSOsvVPbaxHQf3dELf00n06"
    "tioavssTdQzEZgXYOh1AyqtSSJkuA/LZ74qwNsLxvLHDNo5qkOUBp2PmR09wTy0NEPqtNh1I"
    "F9L9+tzKf0udyUrm21XAzuwWOrpKx4O+nYr9yXY8Z3qO44zoBPEg8f8IMUYqcW2ZLTuTDUny"
    "jRQANw0/A94e4k/sKFlyDdlkZccKz8lGBsoXDeWZCdL60aX/lnLF2EiWEB/LwWHsx8fboeil"
    "PhjGEAAsoZW4rzP/ixtE7FoIi7lF8crGrgHScXHw7Ng3cBuBP7iDyIzeS6wGkPfFJQ7IpySB"
    "Ow/ivD8e/VGschiNNrNwUAM3YLxhmYa46V49hAeE/clS57ZfF4b1mbMpbaOExz7ARDMjHsKj"
    "DLxfJw3nSf7CHcmtdQ/Ni0PByi1SjW4QZeOvhLOyz/Mfc3OVwO5Mz8w8yK0vE7XgG1IpfEx0"
    "XzG76fLBPHX1fUUKRMh6bMLxJBRI0xEOK+9OCB1fFTLsv3MHYwHbry3yckiRVi6gGbOliPQa"
    "/87U1o8ngJHvjJmFKH0L4G8Jsu06Xeisp9s2p0ZobHexhrxAjNJ6xns2ulBfmT8MAbYNResb"
    "0t0Y0GizovbfuaODw3ai5kurDC/7QukiTdL+smg7wNfx8foX5wTQsaFvv+spZ1ICbSDDJKw1"
    "vywglEWDePwoP6o6E7ZnwFXrtYUXRrw0npnqwCAJ6OAWCPO137nDRTSMgQYhlrNxPxBs5JgH"
    "kPVBrvUOiJ8WWXa07nM6bVIeqihHB/+wWt952kdxhCt3MBEpTnr79ufhdYhZ9C3FJpWnj+jA"
    "IqJZEAk9J0mG/c4dgzjwt+gYe7uZbYgbTC9+hLmPGYPCIf6Px/v/LuNC767g2NHMQT2onvjn"
    "vLFZmcsMfHoE9PA6ZokbI8Ksf29ouTJYaoH4x7xJfDHW2GkzE0EofPmndhBmMcUDE6XWDU5L"
    "gIiaTMDNqxraLp/r0+s/0nLZXcNxQlOgXiNvFvL+LmyAJQR6AuLigYsNr8T3WdLjfmmI5JSD"
    "UK4AiHEQHut1JjcohAUc+VU7QgKhkmwgekbreNeOBrOBootNm/fL8gssfFBmDFb11qD2a4KR"
    "J5tOuvRizJQvoSRFTpW5qgpIA0HXad77UQs9gnUtHy9U5lFBRDmTo6jSZ9XsV+3w4CVZWu+u"
    "XICf2mHUpaTjNZBPrWpyqA/L0fGp+HUiOePWQth6cIPMrNZ2bKWtbD0LgxCPHhXJuFns6Md5"
    "nxXcvjV0A/2FptIRC9dtRYOBep4r/Kod700bsb6LPqhMv2vHPYtycgw0jQP57Oqn/BQvZ/0P"
    "mkXAchL+wH5QhhimbkLfW6CuXGdbFXuhq4eSZxqj41nbA3ZSn1cnG4aHCntGZbBtMe/eAYx7"
    "CwLdd74HA0z/1TuQHTeoJiSR5/54+mPa+MPQMJ8LgY6ebt32ifPtJhH62nXFQDVzQ+gUQ9Wx"
    "bZzxHzhIGIPjZWbx77nGdAySzjxQSlr/9I6wQIOP75D5yNz/6B2huxY0nUt8ro8jYA4XfRdh"
    "n2sRUk7i/6Anl35JVSHCa/JXAYCBTIybWtf1RJgETkuVwaUF98yhVeMGDKOcz8T3/d07tJpn"
    "zBLvTH5hKF3lr94hQmp26CjRZvLH9R+jv7n0XLfzQuUFfZJBdUj3UqGkoBEGzgIA1Wfr95ju"
    "Gk0f7guoPDeHDE+LtzrI7cpb9202de129o7dxzszjua1Pcj87ncd6ad3jG4e6Puv//j6j5cE"
    "pKQzcEv+zk2ipLalg6ire/MuAHQLriKhA/NudJoaPxPg641kafGwYsxDNrPzPbDKRQmzGaAe"
    "rR7VDoUsgKUb0a5PyAqynPUwuWj+dofLRxePkjsePbrv9U1WJaUT9vebyqqIcvynAMDkwjSd"
    "SBgNHThy5NnUBkvsjYDJeLrtQRz0OsoyDdoRZcAuqawB192fME48Z53r5IP4mSeIpsruzTaj"
    "6YclwcNHzDHW1rdtfe6hXmqubu3SvdNT/TAMQ3oBi8ftTFiGM/2cyFWD9oRNO14F4v5eFX5Y"
    "Y7C9joABYQEa6HYDR0gFdSLh5w0xivNrTtdL/VSCPyyI2edygz3u3I6GWH02Q0IQVzbbuwCQ"
    "Rt8XqFzuM5ZtezQhXTn/4but19xKNG7pFNgTNUrTc4R3gtxeDKpEn/doqA+CjfSMevaCu7aj"
    "3/04/5XgHFDrlF2Xep0X8PO6MbYbeKXifhcA/LVKOCNjviWBz74TrrdjRntk85cb3d8DHbq9"
    "bx33iEB3xTCJUXNQr+O5EppfFcyBziA/CDN5QjLEkHt8vv8FNbOnuId9yz54e3EoYb+y29GC"
    "YaE/BYCO0P5RkyXyp8xswaz2NPSCpM+CeG1XSdeGgEftr6ZD6BrS9OwxEuoSkgjbEmvXUdb9"
    "jDNpSmgb3CzH/4D64/qJGku6mlKI98XE8KIVxMLI9shPAWD6yOeFyrK7ho88IfONWxCeuE53"
    "2fS2YcTc+LaiWoCOwHiJXFJ0dpoB0l5aSu3dYVwoAcoeyFqZUEWWj+v/7iAxipreowWhaI7g"
    "953seQYw91MAkEwhyHkOzVEDUA/MnhDtI1JA07EmNK9hnzkQAicyyQGexIvgtkkVrEXHOFjJ"
    "+Ely1cQKNKgTlip5nv1iH89/i8u80xovI4kNeLDd0dw7xjJSfhcAqosB9eIZ1uFPN8/tomjv"
    "k9WYVY7zXginawT0DbuapeOnKOS+oCyliJ8yGIf81ynPQwf3OijZkDuXHFEzPr3+NOEp+iWI"
    "+dRiNu4XQjgB/VygFB+zAHC19ZrJ7KtlPOq67VPpuRCQgtjs2ivTanPwxHCMhLgI3yU8Jhl0"
    "ezM/jKMIrHxOBilwNxFimdQCf+7j6T/UYaRp5EQTtVdsCH+SFgGhvfCIWJefAsBa2j47dfid"
    "KaRrbwMpI1fhyM1Tmm6uY1K9ePSUe1vAc1h2MaSsOTWJEV+sGqwwS+kY9cEYihG21Zk32j6e"
    "AFRwoTWHi7jZtKRsGjOlU/wi2J3qTO69iFiQ6oXnnatb4TVt9qH4Dgy6v1EAPSJ1ffaRxnDP"
    "mCp4jWL21Ym67uOX4yNpTSuz+UC7WiGQCf63z65+auDSWZTdrBUYkaG00iQePzWKlaBtBnTq"
    "dYhdIIcljkCO992FOg40aDjbg7iYobt0dewXM8A7+grOkU+kMUEvcou/BL6ZBQobxhHPUio1"
    "wMf7/8vsadwmaiMEWR4yOrokWggoYa1k5kDfPid6Cp4UBoTXTBCsr7Os2wIX64e2qb02WpDR"
    "wDh8YBvGNt0iAuWMWAEx31+AD3oFJxAN7kYtqfe70Y/7P7D6WF4C8gtBOj8xCKIHO9jMaC9L"
    "GJ5WQif1Bwz8dk9uEh8ZzwRGU/KCvMkM9QbGpOqw78zeUXs9a2g3mcAXTeWvwHdYUflw/Fx2"
    "782Tzk8v/7Yuxfba8bkK9I1OM7fNSEtS8MlsikuWIptxHQ/ylB6JXlfcBLNogbwxd3T5HuOg"
    "C2hABwKnrNEz8GUSHzb+TnyWkhe2wamLSTt57o/zPx8DOHRbBoNb6SGRC/qltSQsH86uTK23"
    "ZZYijwV6puUlSd6GQepr3MwXEVLkbCEzdfo44NqBeRPf6z8TX55Xxem9KYNBYkPS9en1T/kh"
    "cnq/hGGipDVTsc1u1pejs4gRI8IUPP00M3mP3DYiqhWg0lL96tH034NDgYJRBOW/Jj64W4+8"
    "IwpCAEjNx73fe3ahZeAF12tPw9dUyWxxKI9VSAPwzbVojw8Mu92UOBC6LEB0sLX2yMPVgkzb"
    "e3AItBmV/B+JL9gqy0wijRRkX3kMH+9/n2ssNO4LR8yW/dFiRD4swc8ub2sSIv1EO4Z8N5Zb"
    "LhUctUTWQ+0XQZyfEeQjiWnH5uls//yvic+foUnWrNAW8gji894fRL9xvV0r3hhlRQmV8pZf"
    "qy0toJmDpgvasGOpHJuz6OeAXvi/pUz0EphxsTF+EesQQ5DfQ5P/lPieQ5M5oY4IZ06NEeTz"
    "/f/7GpP1SMgEOEIWa2jq56tKwY4jWqQtYPpWgW+nmU3LYSA5chgRFyQAE+7VuhQDWi28aPNr"
    "aPIfCh8/Q5Mktwn7XpbxdMSP9785ZCiROBZQ3YVd2raao9d3WxKiAXdsGOnPO7WMZJXUbpfX"
    "hvRvzkur6I1k+QxIGqbehChE+q+Fr5+hSW78ScwgTe/j/F8oAPmBvA4Z8Bqckhju8DUpNhJI"
    "L/b1zFnNMYe4ILFRUuaMax8sbsvW+1hIva0GyonwDpGDyss/FD7/GJpkZpMEAecmNrN//Py9"
    "XkV/FUqWbYsSFKrpdN7Ie6VDl7WbvcxDrAJjYL3u2TDKhXYeNR3Dwng85IPzXDlZArfd/2Ph"
    "+9fQ5H0x2jA2Ite0IdaP85/rOepkbDonlgz7MUgiwTxITrYCJl0LxDXP9o82tjnHIRZJ7TE7"
    "IpDJHvjuWXhBz9dLLZd59X9tfGh/H5oMZBwNoiJd8M/X/9vruQhVuS5ha6tnYmJ3MjSsjab9"
    "mIPAai25IFEOqszCAE9kli3WBNbBOk6KFAlkR6eXy6VN2f6l8eX496FJCVb4Rz2zV/h/IQFy"
    "Numbd9FIM/OxGLsW+9JwIvEd19uLFwwBuaGCoyNnNip4pTkf8K6E72t7SJCuPFeQqPYI7dxC"
    "FlHfjU/nvw9NVgQR+YV7S2j1n148zEZ/FYlXDR085LVMwIbH/Tp3JHywb1mAnC1RXTwTyqvN"
    "2iHhIeWeufvwRs8ecUAQfTNmoVL4JR27mI1vFcS/D02Oo9AGcq9E9fLx/g8ry0587FnNWfyZ"
    "jjb9ahuXcgMx0TEVazT4+mknWMkZ/GaDXDrcZa7evPcg3H65UDma5dIx7d+Nj7MK9h+GJjeO"
    "OFGhYXBl9cfx74bo9og1IDlvc6ZN2nmXCfVLBC3R23WKpHUWOebcB0JkeDdIh1aZvtbYJqZf"
    "D6ivnSFD8qNsARhnTA4g/zA0ibF/t3lT9wKlfXz+cdmz3mvQ8OwB2frMYq5zOgFmuicv0PyC"
    "wA4d47yzQCH+XSW5g9x6I9c9xEqkc8dgM5d/VyBlejyNUElH8g9Dk4Ku+zCoQOg07cf7vwsD"
    "1d4e+zW4AjVntZV4/2OO7VS/R/Tc+1UZ9COvUtQbQ0PGP3RkeMcc9Ib4TGCMxoE4p/Xr6WRn"
    "c1TiPw9NNn0sDAJfnZqTIB+WXIJr2awE3viebHTOhGyvc6CLOm0iMtfjNbdiAWVcXQhc8gzL"
    "m9zke3hh30xvuYtR039sUHdLN43s6T8PTe6liQBeYSzVH1/+bGIo1MAxhz/xv+uDBu3zDs8z"
    "kx2E3YxeN6Lb9jrwEIXL3oPDw166dXOsz5pxQrk4KsGN6GiAR3iMH7BZ/g9Dk201AoNNfu17"
    "Ux9nwDlu6JFSWJYdQ31b+auLF59oB0/OdEOblzEjVzPoByqa+zo7vSZfGIdHFNvbgrQmnEh8"
    "id3Q4MHoNYJMkYn/PDTJg+/yXGIFpvvH+7+GEZdEP11mTXtWNiqCU+Q8h5vZ22WZjTAsoCGr"
    "2A1BtMvYvrzn9oXkofaMS7gIn22knG2dwcbfjcNyi529T/dvQ5OtpJr8vDKJCggf93/W4SOD"
    "w3AnJLRGkMu/QCHSezCeF1aEEaZZV6nYwm9lrSypiieqi0gnur/3YOdy/THO4troFYMjms2/"
    "D01SU5Ya3RATWbqP33+SWkId0GjEfJZ4srdI80ANNttZemlXH2yEd1ETwQwRHOF9gnlxDxdz"
    "4K3ssyFgq7Mffnkjoi1PGN0L1ZGq9rehSaJYlfeQbdbLERR/vP4H8ajMec/xgdH1n3zv/Cow"
    "b0CigRtd25OJXihgUA8RynHtq8KDdratZWa3AenPdu4nmk9BPUKA+x6Mg92CcOTvQ5NKIwq8"
    "qBAM1p6ej6f/cZXmNbENUtHD7he6gOuBd1Ym7YUpDNSpg9luQHBv743nsl3dzHszrHa2Ogv6"
    "DhjH+rWG3sNZkejNZiphV+/SX4cmJwpKazBupYmir0S4eOiP+38LlFwvSJPczMlEDOF1A85x"
    "D1qWXNqMRyvllbVYC3/sWqVUPnonETf5UYeBcRGbhLmOvrnJjO0CI0viUi7yL0OTuwdW1txn"
    "x1HXyKyo5enj8x9cC+IQ7GC4tz9k3NsXMXmzlOV1Tds2xrU4WlhdOMP4XnCFqndR6xZFvucN"
    "JgjvjIetMRZmchNSmgPBS2n78efQJBBHpBbOE9Pw1N2cnY/bxwHQlRgejK/waDMngcCuwviU"
    "t5MGx3u8HBQBsZoeHjs71n5GoPZL7jM30GuaFJbMdTwIcPa1ZMqO5eiIK0OofxmapAiZDI1S"
    "4Q+R9016ucaP5783GyluANKACKnmBPbUIGxFAw5HHRt5zWy9hzoSzJH/SY3e7ZJvH7FC7DxB"
    "XI6Mmlw2j2Tw6P1GpuBxH+DPocmFUYlb4rUxPGuo7t1Owz7e/5dTJXzrgs7Qle9zAVR1xmxl"
    "wfWSYppBfUG46+btFp7NtP4x4/0bMMBBex/JS/mTypgbFNO6vHRq0Qfyx9BkFkxJPXKeCREP"
    "olBSZ/P7x/NfTGK4UrOj6Q3FnusQbD+r4pCUnikhsNZbq4lGwuYIb9bnC3dpJgJrXpRDVih0"
    "QHD8VzLT97IO83to0niBSJdHUm6yBM2JjGURBENi+ngF1ImwgarpNkfBs6n3HZGsjVGF1mQy"
    "N1zM2KtknFORG8k9XLtGAqdmKrww6ZEdA9ujANwOT1ADkPrHNShyhFrfmRN4UZEQWhY+CKV+"
    "R6BBZR5OLfXj+f9qWfTcN5fSvm47+m4/07kiULeveNJ9Foe3lRoWEB0v4E7k9hgA3lc63Yom"
    "tJfXvobZOngiDOqtpdGDEDuGxFLnFO2OlLkXDIGuY+SbhdGZ9bHx3BX9/P0XRWxtR8KnYT2P"
    "CxdoCPIWwqhCR1/mdYWz11luWuyrrUZZcyD0Vem1IhV6TRsmyzrL3UduuAHPde0u9URYiRqD"
    "yTVYbhQcmsGh9gKbO959ttSrJVhPP71+Mib53dgc7rgHRnJqaqIRGKIdhTiImwt5QcrG5Bcq"
    "sVcQCRGhsxOJgKnSEEmQ0hGY9wSTOS+5p3WCYin1gVqzbBg66wxz4bwOuSA4sgg1wMBK9Zo+"
    "fv9ptIGcgZDQ85hJPJBrne0OwrYNiNmk416iU9d4mluL6Aey1nMOgK1HRBe44RbA4yiGACuJ"
    "lyJFo7mzSG7WhkFfm+FcRrALWvm92Rkl0swbi5LE0j/e/zRgtQSsrHed1x5fe9k3oRwcErkQ"
    "IvTdMKtZ7QbxrkCTZn2YpbbJ/+fFUEVqr23I2nY671HIHh2IvwTv0t5yTr6vW3fM9J164Cr2"
    "sYo1HAiLYz+iah+f/+UYlKyUZp03tbWXP0tf0RpQndEnLCBzWihvVA18kerDk1wtJerolJL7"
    "aISS7HmDwfjF88pcCWNLLxcJy6dZR9S72pD+ho0S0XomYyIMKscoLN/Rf9z/t3ntRZ9xKJp5"
    "B5hb9byyHHFg5WGgN1jEvN3gfhD/wf6kvlKupdAv5sl7aJJohfHMIqZn+MMaET13CJiO992g"
    "+9WXiIqEP/rT6f/MtpF1Ek4daHvcZxcP8/o/dHGqnoht7SzlonWiW/dZwvPab3T/BqEr9IAU"
    "IatoZtrnLjJd7N25P4cmlZx3QeFSiLS+RsPEvuu2vhFVZa2Cqwcl/Z1kz8tsAhuzafiBi9r+"
    "cf6XTXMm5zaZWJt3Fi0mzh4WWe2+hTMopa2ZRzmRrHtj14HM1qzHvw9N5t07o6Kt6Rx23vD6"
    "gG6BIpfOCAHtYrUduSkEvTyD177N3PGHZV/wMbYVHfyccOjo9+d996sxMfTdRiOR31lYg4Fw"
    "FaRxFBpdl9xzjn8fmixbwiUqJhyhBrFAgx1EvGbzw9K5QYfZmWZzlAy9yyyog94+v/4zWc8c"
    "1JUXCDvnOiNoRUys151bAVJPZIvKEV5H6ZpBjcupZt9+WSH9y9DkReXqGPEIbhe3DvT8MK9+"
    "xeAvq0EO3fKBCpZL5W33ggGxED5e/91XWaJxhiK1ARITpeI8GAjRhkaKss7rKmMHub06Gnjb"
    "d4R8pM2ed62XJf1laFJnsOXY+gHm3OZkvznntPzMlarLw3aeM8B2DURnmY1o5z4+P//yM+mJ"
    "aJ9ZRGuQZ0PjKAPKuRDCg6rUlY3011PJAbeGrNScfOgNETJRwfw5NKko8b0/T0cUlVEzNIUN"
    "ZutjY7O2UG9wA1SAWWGDllcooz4fx/9ArXTjWDSIYPBMR6bZnnCVCIvJhONh7+OaxbBsHlyk"
    "WzmCY/syNvPiVQ5/DE02Ziy6ivK8ywAnmxekEYUGnkPQ1vE0+Gk8RPduBLLvoSP4ePyX0LMN"
    "SHo1574PW6oKsl+pz8G36Bu0UXScwW2Jdk7LQ1/M8WCgh3jo0fzifg1NYggNcwAW1xRQRXi7"
    "hsfYhzviwPdjV8EXjCpuXAKY1j+Z/4/Xv3aDOk8I9bEzQGa+H4PC0lLPJsZl2/L18x0V78dt"
    "BZZbbdmcQweEh+o1Zhco/AxN1uTW2U5pA7+OWVjQeNCoE6Xm1T2nNAp5xEgYT5E85J4wfJqP"
    "538cEzP0pcwQCMxb//ZCCTp/ZDGRIlrZTyQrS3j3acySPe9zmOVKuP6A1GemiMgMBX7faVtS"
    "eieGGLyaB8ZHFZ4jr3aRl33aPqU/V35wH69zz6A/nv9rs95B99dLw3LFtcTFzmtAlknwfD5e"
    "ePBzuD/9XNXwYCxEG+jk9cySAamMsI77Na8H6Z1XAxeP2/zJXqMT6PjndwuARNMZtU0HiOEW"
    "+FhmXzg8JXweABM4X+yZiXASUPMxhoXj7oRX/sBsbd+DmJOKZj80nv28uzq98syBD5Nfo9SU"
    "diD7jx37TeA7a546cM3Wf7IfDuIcjV/W+eFzatiOcXddJEaHo30c/6IVu3mrDdfX+yxiGCfV"
    "6LBOh87+PdRvufbW9NQwLAr1qMf/urvifpbGTYseg8T7ClmVUrSJpTTiNishj5R9QH51h2qw"
    "Y3SdQ9T64PVQLsVZKP14/9eOj6C913q1PzcSMMZXWEbco75vGwOMG723r4szeg6LgYqAMAh/"
    "sBauEMFjOKhSo+pHsaJnH5sw4PYTDAKmVJdV6xr48oS9uwSLnXetIi80s97Wj4/3v77uQ75R"
    "YFsFe0+zkwS6Y8hur12VA7YrlXvbe63nvN7VzgtOESGBM5WBPK7ex1btgux5eOksIUMK5pli"
    "si6g6ghsZtbX5cH4Jw6E0sFcINefzs/t4+tndSwQzry3uJp3LS8W9N8z26X5uvHtTrDt4lgo"
    "m2MNg47T4m/1TRFE8JFzyhmiYbcj/CMwe2MNwcjA8CW1dURXQ0IBE6VagEHpzVo2uyzYj+f7"
    "eP0LKFolh7G12Od3gNHA4YpIYgZoVGIy+f48JPfGKmPAvOYIbmv3s5Rf99eQlfCr0Pe/I3tE"
    "K0IQPJkh4sf8Uy+8Z/8Dw49g+DmUrS5eB12fj8OfmcZD7cwrPpnsM++DK5UF/TXG612kBnGd"
    "h4TEcKZqJwpyrzm1vEZEyKwpfjoM4+gTup+XOUdt3OyTeDKSpfktP3MGlnJhRyJ5dlWzgXBh"
    "O1IPDwKr5+P498SDnBcgzEGfXCYX+rmTCv8/jSPEB+xuCdvtMNplZY29tJNkfm+SceW2ra8h"
    "ACHHslBeSCk+vm+168iRLq7EvAiR1LY9SHm7GTe0U7QtTQK9CuE/3v/0OHmjY7bOEZnfp3ET"
    "hHzcIwjeNSL5MtCRC4dstW0jl/1VidHKDrvs/WX8zqTOVobOyGIXTZAUg6TNmAX3akHMYzcG"
    "vlofCuRdPgs0vWdi9grEFf3x9XMJMldScxVLZwPtNt4I5ucNJ3M4cR8bevFUVFuUUptbd8QA"
    "zSlJi5c5+DV4pY7cV2r92g0jlCFuTit6UJLE2pQT4gnBSxBn4rLB3lRFjCwHwgHB+cfrP7Ol"
    "e+leUn+oRN2lPbQEUqV1XnrDrmOvkqezzAelJkQOvASJJ2k3NPhTFctKvRzflI/tJkil5lWp"
    "G0fguxxbEfuC4WNyCMPNpoGKPPqSi6Ee179+Hv6JNH3ahRie7WiisM47r/zybHBBWvC0JZJY"
    "1FoWO3SuUT+EE7H39x0OnvN5me9rMSvGs3U2wh1bq6nM1uiGDOFE9ZljNL/GnNrz0N0qZISV"
    "QiMhfd7/ZT7Hc2FtaKG5/+pHM2Ne5x7mlzh1OfO8tZUb4riI34LPVel5h4dCO2YLIlmQaT3W"
    "RKcLPcriHILBNJHtiiahjpLe13y+Q/2T0jO7xPeaZ13Yfvz+m1dnagZoU0lYVQ6TkSIxQTVG"
    "Hn9yNAbXEnv84dzrQeSX6Wxqn3e4VPDO4ZbddDY8He8vTsGgII1c+6T186tSpXTH+w6YYXwM"
    "xmmozM0+iVQumldvPj7/eIyVz6+8WbzmyHvnt7cAbSwHSrJ7Z2d9yXZ+KepdDxfR5nMhP3f4"
    "6PdYm4mB5uiYHkeXRrClbCE3joZVnNZ8Q27hFmbvs4U6LkBtcSWuweiHlLF/3P/TUgYXdT8H"
    "LpaPOq/oYULrvNa6zMwPRSNHHINnJ3lYq0Tl/3WHU1e65JnHikQpjJgyMdfRtRmJVrWIYWdX"
    "rOBQjrOycY2956vPyJLPCwPNFnOUHz9/wraVQOVnIimq7arnqXNc1lTy4vR73gHqq2YzZ/eJ"
    "bwLR/s8dXhB3Ol7rvCIAld17uRiqZCOzFRghz4Z04H2pLG7GeVdGS3YIj8KEWJQSNJaDfDz7"
    "jUIrBKDorsI4iGk9jy07tAizWAk1HGw9L3hs6vOOd5WW5fcdbrNd7CAKGeArU9vTvCx71Z4A"
    "ry/QlOJWAKH7uys8PA3YzAikrsBvIB6f4t7n6NSHZU5w+V5P//4WvNn5jk92C3FStiCjE3dI"
    "AUYz+92B3z1v/Y87/GB+a5JSzwN3Q9/P7bKUdcKm4xlroWpFmBN8+4lxz6mO1BQEgktWLM8L"
    "4M8qP97//nhr4dx9UZB4wVW56RMGnC9N2/zeA8TC4YE9nQuk1bBw/b7K5j3nipAIHs5eePpC"
    "FsuP9xfe2kt4q6fTQPBbkPLOSZm+1FlCXRZUqqbinpAHmY/n//rRS3EFyS4C4b2AUNbbdxv/"
    "vMPTQUdc9JpXws+LgdjiOfnjDs8yUx6zl+VBXOiTWVyc33k9x6jwR2r3vszpx/XVosJN7kAa"
    "4ox01IK2hHYDRH++/IMOes4rstnMQg7Euly3n6z8vMPVrIX32es2y9trmTZM/rjKptpS319y"
    "/W6dbHxVQc+vEDwRCqK5y3ymsiGCuDu6EsE4mV8x3Gfpc96N+cZDn4f/v+QgCz7qVkKJfuYs"
    "trmuGaDLmF//JmaZ5NVqcPEvV9nUjcp3YQD5TyC8mrBIDBIzydv7/r4BSWCYyPJ12PkVu/W4"
    "MerNpMn7twjIz/f/f+UrX/nKV77yla985Stf+cpXvvKVr3zlK1/5yle+8pWvfOUrX/nKV77y"
    "la985Stf+cpXvvKVr3zlK1/5yle+8pWvfOUrX/nKV77yla985Stf+cpXvvKVr3zlK1/5yle+"
    "8pWvfOUrX/nKV77yla985Stf+cpXvvKVr3zlK1/5yle+8pWvfOUrX/nKV77yla985Stf+cpX"
    "vvKVr3zlK1/5yle+8pWvfOUrX/nKV77yFYD/B92aGZl3Kab3AAAyGGlUWHRYTUw6Y29tLmFk"
    "b2JlLnhtcAAAAAAAPD94cGFja2V0IGJlZ2luPSLvu78iIGlkPSJXNU0wTXBDZWhpSHpyZVN6"
    "TlRjemtjOWQiPz4KPHg6eG1wbWV0YSB4bWxuczp4PSJhZG9iZTpuczptZXRhLyIgeDp4bXB0"
    "az0iQWRvYmUgWE1QIENvcmUgNS4zLWMwMTEgNjYuMTQ1NjYxLCAyMDEyLzAyLzA2LTE0OjU2"
    "OjI3ICAgICAgICAiPgogICA8cmRmOlJERiB4bWxuczpyZGY9Imh0dHA6Ly93d3cudzMub3Jn"
    "LzE5OTkvMDIvMjItcmRmLXN5bnRheC1ucyMiPgogICAgICA8cmRmOkRlc2NyaXB0aW9uIHJk"
    "ZjphYm91dD0iIgogICAgICAgICAgICB4bWxuczp4bXA9Imh0dHA6Ly9ucy5hZG9iZS5jb20v"
    "eGFwLzEuMC8iPgogICAgICAgICA8eG1wOkNyZWF0b3JUb29sPkFkb2JlIEZpcmV3b3JrcyBD"
    "UzYgKFdpbmRvd3MpPC94bXA6Q3JlYXRvclRvb2w+CiAgICAgICAgIDx4bXA6Q3JlYXRlRGF0"
    "ZT4yMDE1LTAyLTA5VDIwOjU2OjExWjwveG1wOkNyZWF0ZURhdGU+CiAgICAgICAgIDx4bXA6"
    "TW9kaWZ5RGF0ZT4yMDE1LTAyLTA5VDIxOjAxOjMxWjwveG1wOk1vZGlmeURhdGU+CiAgICAg"
    "IDwvcmRmOkRlc2NyaXB0aW9uPgogICAgICA8cmRmOkRlc2NyaXB0aW9uIHJkZjphYm91dD0i"
    "IgogICAgICAgICAgICB4bWxuczpkYz0iaHR0cDovL3B1cmwub3JnL2RjL2VsZW1lbnRzLzEu"
    "MS8iPgogICAgICAgICA8ZGM6Zm9ybWF0PmltYWdlL3BuZzwvZGM6Zm9ybWF0PgogICAgICA8"
    "L3JkZjpEZXNjcmlwdGlvbj4KICAgPC9yZGY6UkRGPgo8L3g6eG1wbWV0YT4KICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIAogICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgCiAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAKICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgIAogICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgCiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAKICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIAogICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgCiAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAKICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgIAogICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgCiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAKICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgIAogICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "CiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAKICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIAogICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgCiAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAKICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgIAogICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgCiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAKICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIAog"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgCiAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAKICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIAogICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgCiAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAKICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgIAogICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgCiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAKICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIAogICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgCiAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAKICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgIAogICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgCiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAKICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "IAogICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgCiAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAKICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIAogICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgCiAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAKICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgIAogICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgCiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAK"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIAogICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgCiAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAKICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgIAogICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgCiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAKICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgIAogICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgCiAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAKICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIAogICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgCiAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAKICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgIAogICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgCiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAKICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIAogICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgCiAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAKICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgIAogICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgCiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAKICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgIAogICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "CiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAKICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIAogICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgCiAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAKICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgIAogICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgCiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAKICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIAog"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgCiAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAKICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIAogICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgCiAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAKICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgIAogICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgCiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAKICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIAogICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgCiAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAKICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgIAogICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgCiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAKICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "IAogICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgCiAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAKICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIAogICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgCiAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAKICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgIAogICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgCiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAK"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIAogICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgCiAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAKICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgIAogICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgCiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAKICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgIAogICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgCiAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAKICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIAogICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgCiAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAKICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgIAogICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgCiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAKICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIAogICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgCiAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAKICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgIAogICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAg"
    "ICAgICAgICAgICAgICAgICAgICAgCiAgICAgICAgICAgICAgICAgICAgICAgICAgICAKPD94"
    "cGFja2V0IGVuZD0idyI/PneDgskAABWbSURBVHic5Zt5kFXXfec/d3v70q+X1zu9IQFqgRBb"
    "YAwIsGgT2VhO7EmN4pp4ppxRVI5NylQyVXGVEk95pjJT5aIqi7VEcUpRFMmJx2gkGUtB0qBI"
    "gBAgGgkQvTfdTS/06377ct99796TP95ilm7oplspV+Vb9f55955zfr/fued3fud3vj9JCMFn"
    "iS7tYAuwE+gEVgH3ArVA1U2vzgLXgD6gF7gEvHc0d2jks5RPWm4DdGkHJQoKPwY8DHQssctB"
    "4G3gZQoGWVaBl80AXdrBauD3gf8CtJb+FwgEJjExRUxMEhHjxMU0OdIICmNLSGi48ElBAlIj"
    "fqkev1SHhIKEdP0wV4DngR8dzR2aWQ65l2yALu1gHfCHwBOAG8Aij4FOVIwTFmMkRQhdJMmT"
    "xcKkZJYbBEECJGQUVOw4JA8eqYZKqZkKqREbDmTU0usp4Bngh0dzh6aWIv9dG6BLO6gC3wJ+"
    "APgK82wQE5OMWR8TEVcxyQESFJW9Wel5hSrPeqGtgkZAaqJZfgC/VI+CrfROHHgSeOpo7lD+"
    "bvS4KwN0aQfXAi8A6wFy6ITFKMPWaZJi5hZFVVVFURRM0ySfvys5C8Ii4ZGqaZO3UCmtQMNR"
    "enQe+J2juUMXFt3nYg3QpR38XeAvAYdJnpiYYMA6SUxMzvm+3+9n165dNDQ0MD4+zqlTp5iZ"
    "mcGyrMXKemO/Uj0r5f+AX2pAKSwNHfjO0dyhv1lMPws2QJd2UAGeBb4JkCXFqHWOUau7uK5v"
    "RW1tLd/73vfYvXs3kiQxPT3NJ598wjPPPENfXx9L9T8yCivkB1khb8BecD8APwZ+72ju0NxC"
    "3QT1zq9Al3bQBfwE2A8QFeP0W8eJion5hZNl9u/fz6OPPkogEMA0Tfx+PzabjUceeYTh4WEM"
    "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"
    "iAAAAABJRU5ErkJggg==")

#def initRosetta():
    # Useful function for re-initializing Rosetta and incorporating the new params files
#    goToSandbox("params")
#    paramsFiles = glob.glob("*.params")
#    paramsstr = ""
#    for params in paramsFiles:
#	paramsstr = paramsstr + params.strip() + " "
#    if (len(paramsstr) > 0):
#	paramsstr = "-extra_res_fa " + paramsstr.strip()
#    init(extra_options=paramsstr + " -ignore_unrecognized_res -pack_missing_sidechains false")
#    pose = pose_from_pdb("../data/bigPDB.pdb")
#    goToSandbox()
#    del pose

def defaultPyMOLView(cmd, model=None):
    # Useful function for setting the default view of a model in PyMOL
    if (model):
	#cmd.select("dsele", "model " + model)
	dsele = "model " + model
	# Set to a default view
	cmd.hide("lines", dsele)
	cmd.show(primaryRender[0], dsele)
	cmd.hide(primaryRender[1], dsele)
	#cmd.select("dsele", "metal and model " + model)
	dsele = "metal and model " + model
	cmd.show("spheres", dsele)
	#cmd.select("dsele", "symbol c and model " + model)
	dsele = "symbol c and model " + model
	cmd.color("gray", dsele)
	#cmd.select("dsele", "model " + model)
	dsele = "model " + model
	cmd.set("ribbon_color", "white", dsele)
	cmd.set("cartoon_color", "white", dsele)
	#cmd.select("dsele", "ss s and model " + model)
	dsele = "ss s and model " + model
	cmd.set("ribbon_color", "yellow", dsele)
	cmd.set("cartoon_color", "yellow", dsele)
	#cmd.select("dsele", "ss h and model " + model)
	dsele = "ss h and model " + model
	cmd.set("ribbon_color", "red", dsele)
	cmd.set("cartoon_color", "red", dsele)
	#cmd.select("dsele", "(ss b or ss t) and model " + model)
	dsele = "(ss b or ss t) and model " + model
	cmd.set("ribbon_color", "blue", dsele)
	cmd.set("cartoon_color", "blue", dsele)
	#cmd.select("dsele", "(ss g or ss i) and model " + model)
	dsele = "(ss g or ss i) and model " + model
	cmd.set("ribbon_color", "orange", dsele)
	cmd.set("cartoon_color", "orange", dsele)
	#DNA
	#Adenosine
	logInfo('tools.py 1316')
	dsele = '(resn ADE) and model %s'%(model)
	cmd.set('ribbon_color','green',dsele)
	cmd.set('cartoon_color','green',dsele)
 	cmd.color('green',dsele)
     #Thymine
	dsele = '(resn THY) and model %s'%(model)
	cmd.set('ribbon_color','red',dsele)
	cmd.set('cartoon_color','red',dsele)
 	cmd.color('red',dsele)
     #Cytosine
	dsele = '(resn CYT) and model %s'%(model)
	cmd.set('ribbon_color','blue',dsele)
	cmd.set('cartoon_color','blue',dsele)
 	cmd.color('blue',dsele)
     #Guanine
	dsele = '(resn GUA) and model %s'%(model)
	cmd.set('ribbon_color','gray',dsele)
	cmd.set('cartoon_color','gray',dsele)
	cmd.color('gray',dsele)
	#cmd.delete("dsele")
    else:
	cmd.select("dsele", "all")
	dsele = "all"
	# Set to a default view
	cmd.hide("lines", "dsele")
	cmd.show(primaryRender[0], dsele)
	cmd.hide(primaryRender[1], dsele)
	cmd.select("dsele", "metal")
	dsele = "metal"
	cmd.show("spheres", dsele)
	cmd.select("dsele", "symbol c")
	dsele = "symbol c"
	cmd.color("gray", dsele)
	#cmd.select("dsele", "all")
	#dsele = "all"
	#Had to change this because using 'all' for some reason doesn't play nice with DNA
	cmd.select('dsele','!(resn ADE | resn THY | resn CYT | resn GUA)')
	dsele = '!(resn ADE | resn THY | resn CYT | resn GUA)'
	cmd.set("ribbon_color", "white", dsele)
	cmd.set("cartoon_color", "white", dsele)
	cmd.select("dsele", "ss s")
	dsele = "ss s"
	cmd.set("ribbon_color", "yellow", dsele)
	cmd.set("cartoon_color", "yellow", dsele)
	cmd.select("dsele", "ss h")
	dsele = "ss h"
	cmd.set("ribbon_color", "red", dsele)
	cmd.set("cartoon_color", "red", dsele)
	cmd.select("dsele", "ss b or ss t")
	dsele = "ss b or ss t"
	cmd.set("ribbon_color", "blue", dsele)
	cmd.set("cartoon_color", "blue", dsele)
	cmd.select("dsele", "ss g or ss i")
	dsele = "ss g or ss i"
	cmd.set("ribbon_color", "orange", dsele)
	cmd.set("cartoon_color", "orange", dsele)
	#DNA
	#Adenosine
	logInfo('tools.py 1372')
	cmd.select('dsele','resn ADE')
	cmd.set('ribbon_color','green','dsele')
	cmd.set('cartoon_color','green','dsele')
 	cmd.color('green','dsele')
     #Thymine
	cmd.select('dsele','resn THY')
	cmd.set('ribbon_color','red','dsele')
	cmd.set('cartoon_color','red','dsele')
 	cmd.color('red','dsele')
     #Cytosine
	cmd.select('dsele','resn CYT')
	cmd.set('ribbon_color','blue','dsele')
	cmd.set('cartoon_color','blue','dsele')
 	cmd.color('blue','dsele')
     #Guanine
	cmd.select('dsele','resn GUA')
	cmd.set('ribbon_color','gray','dsele')
	cmd.set('cartoon_color','gray','dsele')
	cmd.color('gray','dsele')
	cmd.deselect()
	#cmd.delete("dsele")

def goToSandbox(extra=""):
    # Easily gets us back to the sandbox location
    homedir = os.path.expanduser("~")
    if (platform.system() == "Windows"):
	if (len(extra) > 0):
	    extra = "\\" + extra
	os.chdir(homedir + "\\InteractiveROSETTA" + extra)
    else:
	if (len(extra) > 0):
	    extra = "/" + extra
	os.chdir(homedir + "/.InteractiveROSETTA" + extra)

def AA3to1(resn):
    indx3 = "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR HOH ADE CYT GUA THY RAD RCY RGU URA ".find(resn)
    if (indx3 < 0):
	return "Z"
    else:
	indx = indx3 / 4
	return "ACDEFGHIKLMNPQRSTVWYOacgtacgu"[indx]

def resizeTextControlForUNIX(widget, leftbound, width):
    # The centering feature doesn't seem to work on UNIX, so in order to center we have to
    # calculate the pixel width of the text on this widget and then center it in the location
    # defined by leftbound and width
    font = widget.GetFont()
    dc = wx.WindowDC(widget)
    dc.SetFont(font)
    (w, h) = dc.GetTextExtent(widget.GetLabel())
    if (w > width):
	w = width
    xpos = int(width / 2) - int(w / 2) + leftbound
    if (xpos < leftbound):
	xpos = leftbound
    widget.SetPosition((int(xpos), widget.GetPosition()[1]))

def startNewLog():
    # This function initializes the logger when InteractiveRosetta is launched and writes some
    # information about the user's system that may be useful
    curpath = os.getcwd()
    goToSandbox()
    f = open("sessionlog", "w")
    f.write("InteractiveROSETTA Session Log\n")
    f.write("\n")
    operatingsys = platform.system()
    if (operatingsys == "Windows"):
	f.write("Operating System: Windows " + platform.win32_ver()[0] + ", Ver: " + platform.win32_ver()[1] + " " + platform.win32_ver()[2] + "\n")
	processor = platform.uname()[3] + ": " + platform.uname()[4] # Windows has some extra data here
    elif (operatingsys == "Linux"):
	f.write("Operating System: Linux - ")
	linuxdistro = ""
	for data in platform.linux_distribution():
	    linuxdistro = linuxdistro + data + " "
	linuxdistro = linuxdistro + " Release: " + platform.release()
	f.write(linuxdistro + "\n")
	processor = platform.processor()
    elif (operatingsys == "Darwin"):
	processor = platform.processor()
    else: #???
	processor = platform.processor()
    f.write("Processor: " + processor + "\n")
    f.write("Architecture: ")
    for data in platform.architecture():
	f.write(data + " ")
    f.write("\n")
    f.write("Python Version: " + platform.python_version() + "\n\n")
    f.close()
    # Go back to where we were
    os.chdir(curpath)

def logInfo(msg, filename=""):
    # This is a useful function that will append messages and optionally the data in uploaded files
    # to a log for easy debugging purposes
    curpath = os.getcwd()
    goToSandbox()
    try:
	f = open("sessionlog", "a")
	f.write(msg + "\n")
	if (len(filename.strip()) > 0):
	    f.write("BEGIN FILE DATA: " + filename.strip() + " ======================================\n")
	    f2 = open(filename.strip(), "r")
	    for aline in f2:
		f.write(aline)
	    f2.close()
	    f.write("END FILE DATA: " + filename.strip() + " ========================================\n")
	f.close()
    except:
	print "The file " + filename.strip() + " could not be opened by the logger!"
    # Go back to where we were
    os.chdir(curpath)

def scale_list(input_list):
    """
    Scales an array from 0 to 255.

    Used for scaling a list of values, particularly for coloring.
    """
    # Code courtesy of RosettaCommons PyMOL_Mover
    mi = min(input_list)
    ma = max(input_list)
    if ma - mi < 1e-100:
	ma += 1e-100

    r = [int((i - mi) * 255. / (ma - mi)) for i in input_list]
    return r

def get_hex_per_residue_message(pose, per_residue_values, autoscale=True):
    """
    Converts a list to a hex per-residue message with appropriate residue
    names.
    """
    # Code courtesy of RosettaCommons PyMOL_Mover
    message = ''
    info = pose.pdb_info()

    # Handy scaling, since we are using hex.
    if autoscale:
	per_residue_values = scale_list(per_residue_values)

    # Check if the PDBInfo exists.
    if info is not None and info.nres() != 0:
	# The PDBInfo is defined.
	for i in xrange(len(per_residue_values)):
	    chain = info.chain(i + 1)[0]
	    res = info.number(i + 1)
	    icode = info.icode(i + 1)
	    message += '%s%4d%c%02x' % (chain, res, icode,
					per_residue_values[i])
    else:
	# PDBInfo is undefined.
	# Rather than seg fault, let's try pose numbering....
	for i in xrange(len(per_residue_values)):
	    chain = ' '
	    res = i + 1
	    message += '%s%4d %02x' % (chain, res, per_residue_values[i])

    return message

def relabelEnergies(pose, allresidue_E, modelname, selectedScoretype, cmd, seqpos="none"):
    # This function displays the residue energy on the CA of the selected residue or all visible residues
    # First label everything, then hide the labels of everything not in selstring
    found = False
    for (scoretypestr, name) in scoretypes.items():
	if (name == selectedScoretype):
	    scoretype = scoretypestr
	    found = True
	    break
    if (found):
	sindx = allresidue_E[0].index(scoretype)
    else:
	sindx = 0 # Default to total_score
    residue_E = []
    for i in range(1, len(allresidue_E)):
	residue_E.append(allresidue_E[i][sindx])
    cmd.delete("labelsele")
    i = 0
    for ch in pose[0]:
	for residue in ch:
	    chain = ch.id
	    res = residue.id[1]
	    if (seqpos != "none" and str(res) != seqpos):
		i = i + 1
		continue
	    Elabel = str(int(residue_E[i] * 100.0) / 100.0)
	    if (chain != "" and chain != " " and chain != "_"):
		try:
		    cmd.select("labelsele", "model " + modelname + " and chain " + str(chain) + " and resi " + str(res) + " and name o")
		    cmd.label("labelsele", Elabel)
		except:
		    # Do nothing on a fail (NCAA that has no CA?)
		    pass
	    else:
		try:
		    cmd.select("labelsele", "model " + modelname + " and resi " + str(res) + " and name o")
		    cmd.label("labelsele", Elabel)
		except:
		    # Do nothing on a fail (NCAA that has no CA?)
		    pass
	    i = i + 1

def recolorEnergies(pose, allresidue_E, modelname, selectedScoretype, cmd):
    # Useful function for recoloring all the residues of a pose in PyMOL by energy
    found = False
    for (scoretypestr, name) in scoretypes.items():
	if (name == selectedScoretype or scoretypestr == selectedScoretype):
	    scoretype = scoretypestr
	    found = True
	    break
    if (found):
	sindx = allresidue_E[0].index(scoretype)
    else:
	sindx = 0 # Default to total_score
    residue_E = []
    for i in range(1, len(allresidue_E)):
	residue_E.append(allresidue_E[i][sindx])
    residue_E = scale_list(residue_E)
    i = 0
    for ch in pose[0]:
	for residue in ch:
	    chain = ch.id
	    res = residue.id[1]
	    r = residue_E[i]
	    b = 255 - r
	    g = 0
	    try:
		if (modelname == "nomodel"):
		    # This is to prevent all those "selector-errors" from showing up in the output by attempting to
		    # select from a non-existent model in PyMOL
		    raise Exception
		if (chain != "" and chain != " " and chain != "_"):
		    cmd.delete("colorsele")
		    cmd.select("colorsele", "model " + modelname + " and chain " + chain + " and resi " + str(res) + " and symbol c")
		    cmd.color("0x%02x%02x%02x" % ((r, g, b)), "colorsele")
		else:
		    cmd.delete("colorsele")
		    cmd.select("colorsele", "model " + modelname + " and resi " + str(res) + " and symbol c")
		    cmd.color("0x%02x%02x%02x" % ((r, g, b)), "colorsele")
	    except:
		# Model not defined yet, but we still want to get the updated b-factors so don't exit
		pass
	    i = i + 1
	    # If total energy, save the color as a b-factor so the sequence window can display this information
	    # at any time
	    try:
		residue["N"].set_bfactor(r / 255.0 * 100.0)
	    except:
		pass
    return pose

def appendScorefxnParamsInfoToFile(filename, weightsfile):
    # This function appends the score function and extra_res params data to an input file
    # This needs to happen before sending data to the server because the server doesn't have
    # the local files, so they need to be encoded
    # PDBs also need to be encoded but they are currently added in as part of the protocol
    # code because some protocols can have multiple PDBs
    f = open(filename, "a")
    f.write("SCOREFXN\t" + weightsfile + "\n")
    f2 = open(weightsfile, "r")
    f.write("BEGIN SCOREFXN DATA\n")
    for aline in f2:
	f.write(aline.strip() + "\n")
    f2.close()
    f.write("END SCOREFXN DATA\n")
    # Now get the params data
    curdir = os.getcwd()
    goToSandbox("params")
    paramsfiles = glob.glob("*.params")
    for paramsfile in paramsfiles:
	# Do not send over nucleotides, they apparently already get imported by default
	f.write("PARAMS\t" + paramsfile.strip() + "\n")
	f.write("BEGIN PARAMS DATA\n")
	f2 = open(paramsfile.strip(), "r")
	for aline in f2:
	    f.write(aline.strip() + "\n")
	f2.close()
	f.write("END PARAMS DATA\n")
    f.close()
    os.chdir(curdir)

def sendToServer(inputfile, remoteServer=None):
    # This function attempts to send the inputfile to a server running PyRosetta using HTTP POST
    home = os.path.expanduser("~")
    serverpath = ""
    if (platform.system() == "Windows"):
	fin = open(home + "/InteractiveROSETTA/seqwindow.cfg", "r")
    else:
	fin = open(home + "/.InteractiveROSETTA/seqwindow.cfg", "r")
    for aline in fin:
	if (aline.startswith("[SERVER]")):
	    serverpath = aline.split("\t")[1].strip()
    fin.close()
    if (inputfile == "testinput"):
	if (platform.system() == "Windows"):
	    f = open(home + "\\InteractiveROSETTA\\testinputtemp", "w")
	else:
	    f = open(home + "/.InteractiveROSETTA/testinputtemp", "w")
	f.write("THIS IS A TEST")
	f.close()
    elif (inputfile.startswith("kill")):
	if (platform.system() == "Windows"):
	    f = open(home + "\\InteractiveROSETTA\\killinputtemp", "w")
	else:
	    f = open(home + "/.InteractiveROSETTA/killinputtemp", "w")
	f.write(inputfile.split("|")[1])
	serverpath = inputfile.split("|")[2]
	inputfile = "killinput"
	f.close()
    if (remoteServer is not None):
	serverpath = remoteServer
    if (platform.system() == "Windows"):
	f = open(home + "\\InteractiveROSETTA\\" + inputfile + "temp", "rb")
    else:
	f = open(home + "/.InteractiveROSETTA/" + inputfile + "temp", "rb")
    # Added support for a server running on the client's machine, using the keyword "localhost"
    if (serverpath.lower().strip().startswith("localhost")):
	serverlocation = serverpath[serverpath.find(":")+1:]
	myID = "0000000000"
	try:
	    if (platform.system() == "Windows"):
		os.rename(home + "\\InteractiveROSETTA\\" + inputfile + "temp", serverlocation + "\\jobfiles\\" + socket.gethostname() + "-" + inputfile + "-" + myID)
	    else:
		os.rename(home + "/.InteractiveROSETTA/" + inputfile + "temp", serverlocation + "/jobfiles/" + socket.gethostname() + "-" + inputfile + "-" + myID)
	    response = "InteractiveROSETTA Upload Successful"
	except:
	    response = "Failed"
    else:
	datagen, headers = poster.encode.multipart_encode({inputfile: f})
	for line in datagen:
	    if (line[0:2] == "--"):
		myID = line[2:len(line.strip())-2]
	try:
	    request = urllib2.Request(serverpath + "/cgi-bin/jobupload.cgi", datagen, headers)
	    # Actually do the request, and get the response
	    response = urllib2.urlopen(request).read().strip()
	    f.close()
	except:
	    f.close()
	    raise Exception("ERROR: Failed to upload protocol inputs to the server")
    if (response == "InteractiveROSETTA Upload Successful"):
	return myID
    else:
	raise Exception("ERROR: Failed to upload protocol inputs to the server")

def queryServerForResults(outputfile):
    # This function looks to see if the server uploaded an outputfile
    # If it did, download it and generate the PDB files that the GUI looks for
    curdir = os.getcwd()
    goToSandbox()
    localfilename = outputfile.split("-")[0]
    try:
	if (serverName[0].lower().strip().startswith("localhost")):
	    serverlocation = serverName[0][serverName[0].find(":")+1:]
	    if (platform.system() == "Windows"):
		f = open(serverlocation + "\\results\\" + outputfile)
	    else:
		f = open(serverlocation + "/results/" + outputfile)
	else:
	    f = urllib2.urlopen(serverName[0] + "/results/" + outputfile)
	f2 = open(localfilename + "temp", "w")
	for aline in f:
	    f2.write(aline.strip() + "\n")
	f2.close()
	f.close()
	# Now generate the PDBs from the data in the file so the GUI doesn't crash when it doesn't see them
	f = open(localfilename + "temp", "r")
	readingData = False
	for aline in f:
	    if (aline[0:6] == "OUTPUT"):
		pdbfile = aline.split("\t")[1].strip()
		pdbfile = pdbfile[pdbfile.find("-")+1:]
		f2 = open(pdbfile, "w")
	    elif (aline[0:14] == "BEGIN PDB DATA"):
		readingData = True
	    elif (aline[0:12] == "END PDB DATA"):
		f2.close()
		readingData = False
	    elif (readingData):
		f2.write(aline.strip() + "\n")
	f.close()
	# Now make it visible to the main GUI
	os.rename(localfilename + "temp", localfilename)
    except:
	# Not there yet, don't do anything
	pass
    # Maybe there is an error report, let's look for it
    try:
	ID = outputfile.split("-")[1]
	if (serverName[0].lower().strip().startswith("localhost")):
	    serverlocation = serverName[0][serverName[0].find(":")+1:]
	    if (platform.system() == "Windows"):
		f = open(serverlocation + "\\results\\errreport-" + ID)
	    else:
		f = open(serverlocation + "/results/errreport" + ID)
	else:
	    f = urllib2.urlopen(serverName[0] + "/results/errreport-" + ID)
	f2 = open("errreporttemp", "w")
	for aline in f:
	    f2.write(aline.strip() + "\n")
	f2.close()
	f.close()
	# Now make it visible to the main GUI
	os.rename("errreporttemp", "errreport")
    except:
	# Not there yet, don't do anything
	pass
    os.chdir(curdir)

def getRecognizedTypes():
    # Returns a list of 3 letter codes that will be recognized by Rosetta
    recognized = ["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN",
	      "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR", "HIE", "HID", "ADE", "CYT",
	      "GUA", "THY", "RAD", "RCY", "RGU", "URA","DA","DT","DC","DG"]
    if (platform.system() == "Windows"):
	# Windows has all the metal ions by default
	recognized.extend(["CA", "FE2", "FE", "K", "MG", "MN", "NA", "ZN"])
    curdir = os.getcwd()
    goToSandbox("params")
    params = glob.glob("*.fa.params")
    for param in params:
	recognized.append(param.split(".fa.params")[0])
    os.chdir(curdir)
    return recognized

def cleanPDB(pdbfile, acceptNCAAs=False):
    # This function will look for and remove duplicate atoms
    # It will permanently modify the PDB that the user loaded
    # This shouldn't cause problems for other programs though
    # I am also going to give blank chain IDs a value because it
    # is causing too much trouble with some of the Rosetta protocols
    # to handle blank chain IDs
    data = []
    taken_nums = {}
    num_mapping = {}
    takenIDs = ""
    # First we have to get a list of all chainIDs that are currently taken
    # We need to do this so that when we are cleaning up the PDB file, we can rename
    # blank chain IDs to something else
    # Sometimes PDB files have multiple chains with no ID and they are distinguished solely
    # by a TER line in between the coordinate secions, but Rosetta needs to have valid IDs
    # in order to do anything with the chains, and BioPython will freak out when it sees multiple
    # residues with blank chainIDs and the same residue positions
    fin = open(pdbfile.strip(), "r")
    data = []
    for aline in fin:
	if ((aline.startswith("ATOM") or aline.startswith("HETATM")) and not aline[21] in takenIDs):
	    takenIDs += aline[21]
	data.append(aline)
    fin.close()
    blankID = ""
    for char in "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789":
	if (char not in takenIDs):
	    blankID = char
	    break
    # Sometimes there are PDBs that have multiple residues with the same residue index
    # BioPython drops these, but it can lead to all kinds of problems later on
    # So I will keep a record of the backbone atoms of the current residue and if we encounter
    # a BB atom with the same residue index, we will assume it's a new residue and renumber the residue
    # Updated to keep the first such residue and drop all others with the same index
    lastBBatoms = []
    altlocs_taken = ""
    #f = open(pdbfile.strip(), "r")
    curr_res = "   0"
    offset = 0
    counter = 0
    while (counter < len(data)):
    #for counter in range(0, len(data)):
	aline = data[counter]
	if (aline.startswith("TER")):
	    takenIDs += blankID
	    blankID = ""
	    for char in "ABCDEFGHIJKLMNOPQRSTUVWXYZ00123456789":
		if (char not in takenIDs):
		    blankID = char
		    break
	if ((aline[0:4] == "ATOM" or aline[0:6] == "HETATM") and not(aline[17:20].strip() in getRecognizedTypes()) and not(acceptNCAAs)):
	    logInfo('Popping line:\n %s'%(aline))
	    offset = offset + 1
	    data.pop(counter)
	    continue
	elif (aline[0:4] == "ATOM" or aline[0:6] == "HETATM" or aline[0:3] == "TER"):
	    try:
		atomno = int(aline[7:11])
		atomno = atomno - offset
		aline = aline[0:7] + ("%4i" % atomno) + aline[11:]
	    except:
		pass
	if ((aline.startswith("ATOM") or aline.startswith("HETATM")) and aline[21] == " "):
	    # Rewrite blank chain IDs
	    aline = aline[0:21] + blankID + aline[22:]
	if ((aline.startswith("ATOM") or aline.startswith("HETATM")) and isAA(aline[17:20])):
	    res = aline[22:27] # Includes residue indx + the optional alternate letter
	    if (res[0:4] != curr_res[0:4]): # New residue indx
		altlocs_taken = res[4] # Reset the taken altlocs
		curr_res = res
		atomtypes = []
		lastBBatoms = []
	    # This is only done if this is a new residue, but not a new residue indx
	    if (aline[22:27] != curr_res or aline[12:16] in lastBBatoms):
		pass
#==============================================================================
# 		curr_res = res
# 		atomtypes = []
# 		lastBBatoms = []
# 		# Assign the altloc to whatever the most recent altloc used was
# 		for char in " ABCDEFGHIJKLMNOPQRSTUVWXYZ":
# 		    if (not(char in altlocs_taken)):
# 			altlocs_taken = char + altlocs_taken
# 			break
#==============================================================================

	    res = res[0:4] + altlocs_taken[0]
	    atomtype = aline[12:16]
	    if (atomtype in [" C  ", " CA ", " O  ", " N  "]):
		lastBBatoms.append(atomtype)
	    if (atomtype in atomtypes):
		# Find a new type for this atom
		stem = atomtype[0:2]
		for i in range(1, 100):
		    if (i < 10):
			newtype = stem + str(i) + " "
		    else:
			newtype = stem + str(i)
		    if (not(newtype in atomtypes)):
			atomtypes.append(newtype)
			break
		aline = aline[0:12] + newtype + aline[16:]
	    else:
		atomtypes.append(atomtype)
	    # Now check for alternate forms of residues (i.e. 30A, 30B, etc.)
	    # Rename these so the each have a unique number, the user can delete extras later
	    chain = aline[21]
	    if (not(chain in taken_nums.keys())):
		taken_nums[chain] = []
	    if (not(res[0:4] in taken_nums[chain])):
		taken_nums[chain].append(res[0:4])
		num_mapping[chain+res] = res[0:4]
	    else:
		try:
		    aline = aline[0:22] + num_mapping[chain+res] + " " + aline[27:]
		except:
		    # Find a new ID
		    lastnum = int(taken_nums[chain][len(taken_nums[chain])-1]) + 1
		    num_mapping[chain+res] = "%4i" % lastnum
		    taken_nums[chain].append("%4i" % lastnum)
		    aline = aline[0:22] + num_mapping[chain+res] + " " + aline[27:]
	    #data.append(aline.strip())
	#else:
	# Change nucleic acid strings to what Rosetta expects
	if (aline[17:20] == " DA"):
	    data[counter] = aline[0:17] + "ADE" + aline[20:]
	    #data.append(aline[0:17] + "ADE" + aline[20:])
	elif (aline[17:20] == " DC"):
	    data[counter] = aline[0:17] + "CYT" + aline[20:]
	    #data.append(aline[0:17] + "CYT" + aline[20:])
	elif (aline[17:20] == " DG"):
	    data[counter] = aline[0:17] + "GUA" + aline[20:]
	    #data.append(aline[0:17] + "GUA" + aline[20:])
	elif (aline[17:20] == " DT"):
	    data[counter] = aline[0:17] + "THY" + aline[20:]
	    #data.append(aline[0:17] + "THY" + aline[20:])
	else:
	    data[counter] = aline
	    #data.append(aline)
	counter += 1
    #f.close()
    f = open(pdbfile.strip(), "w")
    for aline in data:
	f.write(aline)
    f.close()
    return data

def fixPyMOLSave(pdbfile):
    # Read in the data first
    data = []
    f = open(pdbfile, "r")
    for aline in f:
	data.append(aline.strip())
    f.close()
    # Now traverse backwards through the data, and drop any TER lines that are in between two atoms on the same chain
    lastchain = "NA"
    terHit = -1
    for i in range(len(data)-1, -1, -1):
	if ((data[i][0:4] == "ATOM" or data[i][0:6] == "HETATM") and terHit >= 0):
	    if (lastchain == data[i][21]):
		data.pop(terHit)
	    terHit = -1
	if (data[i][0:4] == "ATOM" or data[i][0:6] == "HETATM"):
	    lastchain = data[i][21]
	elif (data[i][0:3] == "TER"):
	    terHit = i
    # Write it out
    f = open(pdbfile, "w")
    for aline in data:
	f.write(aline + "\n")
    f.close()

def fitGridColumn(grid, col, minsize):
    # Useful function for taking a grid and resizing the given column to fit overflowing text
    oldwidth = grid.GetColSize(col)
    font = grid.GetFont()
    dc = wx.WindowDC(grid)
    dc.SetFont(font)
    width = minsize - 10
    (w, h) = dc.GetTextExtent(grid.GetColLabelValue(col))
    if (w > width):
	width = w
    for r in range(0, grid.NumberRows):
	(w, h) = dc.GetTextExtent(grid.GetCellValue(r, col))
	if (w > width):
	    width = w
    if (oldwidth != int(width+10)):
	grid.SetColSize(col, int(width+10))
	grid.Refresh()

def deleteInputFiles():
    # Useful function for deleting all input/temporary/output files from InteractiveROSETTA HOME
    # This is useful if a protocol is canceled because if these files are not removed then the daemon
    # may see them again and start all over
    goToSandbox()
    if (os.path.isfile("minimizeinput")):
	os.remove("minimizeinput")
    if (os.path.isfile("minimizeoutput")):
	os.remove("minimizeoutput")
    if (os.path.isfile("designinput")):
	os.remove("designinput")
    if (os.path.isfile("designoutput")):
	os.remove("designoutput")
    if (os.path.isfile("scoreinput")):
	os.remove("scoreinput")
    if (os.path.isfile("scoreoutput")):
	os.remove("scoreoutput")
    if (os.path.isfile("rotamerinput")):
	os.remove("rotamerinput")
    if (os.path.isfile("rotameroutput")):
	os.remove("rotameroutput")
    if (os.path.isfile("coarsekicinput")):
	os.remove("coarsekicinput")
    if (os.path.isfile("coarsekicoutput")):
	os.remove("coarsekicoutput")
    if (os.path.isfile("kicoutput")):
	os.remove("kicoutput")
    if (os.path.isfile("repackme.pdb")):
	os.remove("repackme.pdb")
    if (os.path.isfile("finekicinput")):
	os.remove("finekicinput")
    if (os.path.isfile("errreport")):
	os.remove("errreport")
    if (os.path.isfile("coarsedockinput")):
	os.remove("coarsedockinput")
    if (os.path.isfile("finedockinput")):
	os.remove("finedockinput")
    if (os.path.isfile("dockoutput")):
	os.remove("dockoutput")
    if (os.path.isfile("dock_progress")):
	os.remove("dock_progress")
    if (os.path.isfile("threadinput")):
	os.remove("threadinput")
    if (os.path.isfile("threadoutput")):
	os.remove("threadoutput")
    tempfiles = glob.glob("*temp")
    for tempfile in tempfiles:
	try:
	    os.remove(tempfile)
	except:
	    pass
    tempfiles = glob.glob("*input")
    for tempfile in tempfiles:
	try:
	    os.remove(tempfile)
	except:
	    pass
    tempfiles = glob.glob("*output")
    for tempfile in tempfiles:
	try:
	    os.remove(tempfile)
	except:
	    pass
    # Remove the sandbox PDBs
    tempfiles = glob.glob("*.pdb")
    for tempfile in tempfiles:
	os.remove(tempfile)

def isAA(resn):
    if (len(resn) == 3 and resn in "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR "):
	return True
    elif (len(resn) == 1 and resn in "ACDEFGHIKLMNPQRSTVWY"):
	return True
    return False