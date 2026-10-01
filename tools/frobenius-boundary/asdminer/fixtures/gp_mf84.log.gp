{S = mfinit([8,4], 0); d = mfdim([8,4], 0); print("DIM=", d); f = mfeigenbasis(S)[1]; print("COEFS=", mfcoefs(f, 3000)); quit; }
