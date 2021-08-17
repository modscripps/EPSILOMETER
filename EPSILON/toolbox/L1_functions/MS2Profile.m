function [Profile] = MS2Profile(MS)
% [Profile] = MS2Profile(MS)
%
% Convert an MS structure (old epsi library) to Profile structure (new epsi
% library)
%
% Nicole Couto | December 2020
% -------------------------------------------------------------------------

Profile = epsiSetup_make_empty_structure;

Profile.dnum = MS.time(:);
Profile.pr = MS.pr(:);
Profile.w = MS.w(:);
Profile.t = MS.t(:);
Profile.s = MS.s(:);
Profile.kvis = MS.kvis(:);
Profile.epsilon = MS.epsilon;
Profile.epsilon_co = MS.epsilon_co;
Profile.chi = MS.chi;
Profile.sh_fc = [MS.kc(:,1).*MS.w(:), MS.kc(:,2).*MS.w(:)];
Profile.tg_fc = [MS.kcfpo7(:,1).*MS.w(:), MS.kcfpo7(:,2).*MS.w(:)];

nFirst = cellfun(@(A) A(1), MS.indscan);
nLast = cellfun(@(A) A(end), MS.indscan);
Profile.ind_range_epsi = [nFirst(:),nLast(:)];

Profile.f = MS.f;
Profile.k = MS.k;
Profile.Pa_g_f.a1 = squeeze(MS.Pf(5,:,:));
Profile.Pa_g_f.a2 = squeeze(MS.Pf(6,:,:));
Profile.Pa_g_f.a3 = squeeze(MS.Pf(7,:,:));
Profile.Ps_volt_f.s1 = squeeze(MS.Pf(3,:,:));
Profile.Ps_volt_f.s2 = squeeze(MS.Pf(4,:,:));
Profile.Ps_shear_k.s1 = squeeze(MS.Pshear_k(:,:,1));
Profile.Ps_shear_k.s2 = squeeze(MS.Pshear_k(:,:,2));
Profile.Ps_shear_co_k.s1 = squeeze(MS.Pshearco_k(:,:,1));
Profile.Ps_shear_co_k.s2 = squeeze(MS.Pshearco_k(:,:,2));
Profile.Pt_volt_f.t1 = squeeze(MS.Pf(1,:,:));
Profile.Pt_volt_f.t2 = squeeze(MS.Pf(2,:,:));
Profile.Pt_Tg_k.t1 = squeeze(MS.PphiT_k(:,:,1));
Profile.Pt_Tg_k.t2 = squeeze(MS.PphiT_k(:,:,2));
Profile.Cs1a3_full = [];
Profile.Cs2a3_full = [];
