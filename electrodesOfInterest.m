% examining the electrodes using an overlay of MNIelectrodes.nii on wc1Onuki, I identified which electrodes where above the central sulcus (i.e. parietal/premotor) and which were below/posterior (visual)
% the A stack refers to left hemisphere, the B stack to right hemisphere. This is shown in color on the pdf Electrodes_parietalvsvisual_ISC.pdf
% based on that I then selected contrasts of electrodes that fall close to ISC Ig0 and are in the parietal vs. visual

parietalA= [58 59;59 60;60 61;61 62;53 54;54 55;55 56;48 49;41 42];
visualA= [51 52; 45 46;46 47;35 36; 30 32;31 32;25 26;26 27];
parietalB=[58 59;59 60;23 24;24 25;53 43];
visualB=[61 62;55 56;49 50;29 30; 30 31;33 34;34 35;35 36; 36 37;38 39;39 40;40 41; 41 42];



