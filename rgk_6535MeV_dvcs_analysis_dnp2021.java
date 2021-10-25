package multi_energy_dvcs;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import javax.swing.JFrame;

import org.jlab.clas.physics.LorentzVector;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.math.F1D;
import org.jlab.groot.ui.TCanvas;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSource;
import org.jlab.utils.groups.IndexedList;
import org.slf4j.Marker;
import org.jlab.io.evio.EvioDataBank;
import org.jlab.io.hipo.HipoDataBank;

public class rgk_6535MeV_dvcs_analysis_dnp2021 {
	
	static HipoDataSource reader = new HipoDataSource();
	static HipoDataSource reader_1gamma_pi0_sim = new HipoDataSource();
	static HipoDataSource reader_pi0_sim = new HipoDataSource();	

	static H2F htheta_vs_p_e = new H2F("theta_vs_p_e", "theta_vs_p_e", 100, 0, 8, 100, 0, 45);
	static H2F htheta_vs_phi_e_rb = new H2F("theta_vs_phi_e_rb", "theta_vs_phi_e_rb", 100, -180, 180, 100, 0, 45);
	static H2F hvz_vs_p_e = new H2F("vz_vs_p_e", "zv_vs_p_e", 100, 0, 8, 100, -20, 20);
	static H2F hvz_vs_theta_e = new H2F("vz_vs_theta_e", "zv_vs_theta_e", 100, 0, 45, 100, -20, 20);
	static H2F htheta_vs_p_proton = new H2F("theta_vs_p_proton", "theta_vs_p_proton", 100, 0, 5, 100, 0, 80);
	static H2F htheta_vs_phi_proton_rb = new H2F("theta_vs_phi_proton_rb", "theta_vs_phi_proton_rb", 100, -180, 180, 100, 0, 80);
	static H2F hvz_vs_p_proton = new H2F("vz_vs_p_proton", "zv_vs_p_proton", 100, 0, 5, 100, -20, 20);
	static H2F hvz_vs_theta_proton = new H2F("vz_vs_theta_proton", "zv_vs_theta_proton", 100, 0, 80, 100, -20, 20);
	static H2F htheta_vs_p_gamma = new H2F("theta_vs_p_gamma", "theta_vs_p_gamma", 100, 0, 8, 100, 0, 45);
	static H2F htheta_vs_phi_gamma_rb = new H2F("theta_vs_phi_gamma_rb", "theta_vs_phi_gamma_rb", 100, -180, 180, 100, 0, 45);
	static H2F hvz_vs_p_gamma = new H2F("vz_vs_p_gamma", "zv_vs_p_gamma", 100, 0, 8, 100, -20, 20);
	static H2F hvz_vs_theta_gamma = new H2F("vz_vs_theta_gamma", "zv_vs_theta_gamma", 100, 0, 45, 100, -20, 20);
	static H1F htheta_cone_gamma_rb = new H1F("theta_cone_gamma_rb", "theta_cone_gamma_rb", 100, 0, 5);
	static H1F hX_proton_m_rb = new H1F("X_proton_m_rb", "X_proton_m_rb", 100, 0, 3);
	static H1F hX_epg_pt_rb = new H1F("X_epg_pt_rb", "X_epg_pt_rb", 100, 0, 0.8);
	static H1F hX_epg_E_rb = new H1F("X_epg_E_rb", "X_epg_E_rb", 100, -2, 2);
	static H2F hQ2_vs_x_B_rb = new H2F("Q2_vs_x_B_rb", "Q2_vs_x_B_rb", 100, 0, 0.7, 100, 0, 7);
	static H2F hnegt_vs_phi_trento = new H2F("-t_vs_phi_trento", "-t_vs_phi_trento", 100, 0, 360, 100, 0, 3);
	
	static H1F htheta_cone_pi0_rb = new H1F("theta_cone_pi0_rb", "theta_cone_pi0_rb", 100, 0, 8);
	static H1F hX_proton_m_pi0_rb = new H1F("X_proton_m_pi0_rb", "X_proton_m_pi0_rb", 100, 0, 3);
	static H1F hX_epgg_pt_rb = new H1F("X_epgg_pt_rb", "X_epgg_pt_rb", 100, 0, 0.8);
	static H1F hX_epgg_E_rb = new H1F("X_epgg_E_rb", "X_epgg_E_rb", 100, -2, 2);
	
	static H1F hphi_CM_poshel = new H1F("phi_CM_poshel", "phi_CM_poshel", 30, 0, 360);
	static H1F hphi_CM_neghel = new H1F("phi_CM_neghel", "phi_CM_neghel", 30, 0, 360);
	
	static H1F hphi_CM_poshel_pi0 = new H1F("phi_CM_poshel_pi0", "phi_CM_poshel_pi0", 30, 0, 360);
	static H1F hphi_CM_neghel_pi0 = new H1F("phi_CM_neghel_pi0", "phi_CM_neghel_pi0", 30, 0, 360);
	static H1F hphi_CM_1gamma_pi0_sim = new H1F("phi_CM_1gamma_pi0_sim", "phi_CM_1gamma_pi0_sim", 30, 0, 360);
	static H1F hphi_CM_2gamma_pi0_sim = new H1F("phi_CM_2gamma_pi0_sim", "phi_CM_2gamma_pi0_sim", 30, 0, 360);
	
	static IndexedList<H1F> histGroups_phi_CM_poshel_bin = new IndexedList<H1F>(1);
	public static void phi_CM_poshel_bin_histos() {
		for(int ibin = 0; ibin  < 8; ibin++) {
			H1F phi_CM_poshel_bin = new H1F("phi_CM_poshel_bin", "phi_CM_poshel_bin", 30, 0, 360);
			histGroups_phi_CM_poshel_bin.add(phi_CM_poshel_bin, ibin);
		}
	}
	static IndexedList<H1F> histGroups_phi_CM_neghel_bin = new IndexedList<H1F>(1);
	public static void phi_CM_neghel_bin_histos() {
		for(int ibin = 0; ibin  < 8; ibin++) {
			H1F phi_CM_neghel_bin = new H1F("phi_CM_neghel_bin", "phi_CM_neghel_bin", 30, 0, 360);
			histGroups_phi_CM_neghel_bin.add(phi_CM_neghel_bin, ibin);
		}
	}
	
	static IndexedList<H1F> histGroups_phi_CM_poshel_tbin = new IndexedList<H1F>(2);
	public static void phi_CM_poshel_tbin_histos() {
		for(int ibin = 0; ibin  < 8; ibin++) {
			for(int itbin = 0; itbin < 8; itbin++)
			{
				H1F phi_CM_poshel_tbin = new H1F("phi_CM_poshel_tbin", "phi_CM_poshel_tbin", 30, 0, 360);
				histGroups_phi_CM_poshel_tbin.add(phi_CM_poshel_tbin, ibin, itbin);
			}
		}
	}
	
	static IndexedList<H1F> histGroups_phi_CM_neghel_tbin = new IndexedList<H1F>(2);
	public static void phi_CM_neghel_tbin_histos() {
		for(int ibin = 0; ibin  < 8; ibin++) {
			for(int itbin = 0; itbin < 8; itbin++)
			{
				H1F phi_CM_neghel_tbin = new H1F("phi_CM_neghel_tbin", "phi_CM_neghel_tbin", 30, 0, 360);
				histGroups_phi_CM_neghel_tbin.add(phi_CM_neghel_tbin, ibin, itbin);
			}
		}
	}
	
	static IndexedList<H1F> histGroups_phi_CM_poshel_pi0_bin = new IndexedList<H1F>(1);
	public static void phi_CM_poshel_pi0_bin_histos() {
		for(int ibin = 0; ibin  < 8; ibin++) {
			H1F phi_CM_poshel_pi0_bin = new H1F("phi_CM_poshel_pi0_bin", "phi_CM_poshel_pi0_bin", 30, 0, 360);
			histGroups_phi_CM_poshel_pi0_bin.add(phi_CM_poshel_pi0_bin, ibin);
		}
	}
	static IndexedList<H1F> histGroups_phi_CM_neghel_pi0_bin = new IndexedList<H1F>(1);
	public static void phi_CM_neghel_pi0_bin_histos() {
		for(int ibin = 0; ibin  < 8; ibin++) {
			H1F phi_CM_neghel_pi0_bin = new H1F("phi_CM_neghel_pi0_bin", "phi_CM_neghel_pi0_bin", 30, 0, 360);
			histGroups_phi_CM_neghel_pi0_bin.add(phi_CM_neghel_pi0_bin, ibin);
		}
	}
	
	static IndexedList<H1F> histGroups_phi_CM_poshel_pi0_tbin = new IndexedList<H1F>(2);
	public static void phi_CM_poshel_pi0_tbin_histos() {
		for(int ibin = 0; ibin  < 8; ibin++) {
			for(int itbin = 0; itbin < 8; itbin++)
			{
				H1F phi_CM_poshel_pi0_tbin = new H1F("phi_CM_poshel_pi0_tbin", "phi_CM_poshel_pi0_tbin", 30, 0, 360);
				histGroups_phi_CM_poshel_pi0_tbin.add(phi_CM_poshel_pi0_tbin, ibin, itbin);
			}
		}
	}
	
	static IndexedList<H1F> histGroups_phi_CM_neghel_pi0_tbin = new IndexedList<H1F>(2);
	public static void phi_CM_neghel_pi0_tbin_histos() {
		for(int ibin = 0; ibin  < 8; ibin++) {
			for(int itbin = 0; itbin < 8; itbin++)
			{
				H1F phi_CM_neghel_pi0_tbin = new H1F("phi_CM_neghel_pi0_tbin", "phi_CM_neghel_pi0_tbin", 30, 0, 360);
				histGroups_phi_CM_neghel_pi0_tbin.add(phi_CM_neghel_pi0_tbin, ibin, itbin);
			}
		}
	}
	
	static IndexedList<H1F> histGroups_phi_CM_1gamma_pi0_sim_bin = new IndexedList<H1F>(1);
	public static void phi_CM_1gamma_pi0_sim_bin_histos() {
		for(int ibin = 0; ibin  < 8; ibin++) {
			H1F phi_CM_1gamma_pi0_sim_bin = new H1F("phi_CM_1gamma_pi0_sim_bin", "phi_CM_1gamma_pi0_sim_bin", 30, 0, 360);
			histGroups_phi_CM_1gamma_pi0_sim_bin.add(phi_CM_1gamma_pi0_sim_bin, ibin);
		}
	}
	static IndexedList<H1F> histGroups_phi_CM_2gamma_pi0_sim_bin = new IndexedList<H1F>(1);
	public static void phi_CM_2gamma_pi0_sim_bin_histos() {
		for(int ibin = 0; ibin  < 8; ibin++) {
			H1F phi_CM_2gamma_pi0_sim_bin = new H1F("phi_CM_2gamma_pi0_sim_bin", "phi_CM_2gamma_pi0_sim_bin", 30, 0, 360);
			histGroups_phi_CM_2gamma_pi0_sim_bin.add(phi_CM_2gamma_pi0_sim_bin, ibin);
		}
	}
	
	static IndexedList<H1F> histGroups_phi_CM_1gamma_pi0_sim_tbin = new IndexedList<H1F>(2);
	public static void phi_CM_1gamma_pi0_sim_tbin_histos() {
		for(int ibin = 0; ibin  < 8; ibin++) {
			for(int itbin = 0; itbin < 8; itbin++)
			{
				H1F phi_CM_1gamma_pi0_sim_tbin = new H1F("phi_CM_1gamma_pi0_sim_tbin", "phi_CM_1gamma_pi0_sim_tbin", 30, 0, 360);
				histGroups_phi_CM_1gamma_pi0_sim_tbin.add(phi_CM_1gamma_pi0_sim_tbin, ibin, itbin);
			}
		}
	}
	
	static IndexedList<H1F> histGroups_phi_CM_2gamma_pi0_sim_tbin = new IndexedList<H1F>(2);
	public static void phi_CM_2gamma_pi0_sim_tbin_histos() {
		for(int ibin = 0; ibin  < 8; ibin++) {
			for(int itbin = 0; itbin < 8; itbin++)
			{
				H1F phi_CM_2gamma_pi0_sim_tbin = new H1F("phi_CM_2gamma_pi0_sim_tbin", "phi_CM_2gamma_pi0_sim_tbin", 30, 0, 360);
				histGroups_phi_CM_2gamma_pi0_sim_tbin.add(phi_CM_2gamma_pi0_sim_tbin, ibin, itbin);
			}
		}
	}
	
	static H1F h2gamma_m = new H1F("2gamma_m", "2gamma_m", 100, 0.05, 0.25);
	static H1F hX_2gamma_m = new H1F("X_2gamma_m", "2X_gamma_m", 100, 0.05, 0.3);
	
	static void processEvent(DataEvent event, int eventCounter, List<Double> phi_val, List<Double> phi_count)
	{
		if(event.hasBank("REC::Particle") && event.hasBank("REC::Track"))
		{
			HipoDataBank rec = (HipoDataBank) event.getBank("REC::Particle");
			int ecount = 0;
			int protoncount = 0;
			int gammacount = 0;
			int e_pindex = -1;
			int proton_pindex = -1;
			int gamma1_pindex = -1;
			int gamma2_pindex = -1;
			LorentzVector p_e = new LorentzVector(0, 0, 0, 0.000511);
			LorentzVector p_proton = new LorentzVector(0, 0, 0, 0.938272);
			LorentzVector p_gamma1 = new LorentzVector(0, 0, 0, 0);
			LorentzVector p_gamma2 = new LorentzVector(0, 0, 0, 0);
			double E = 6.535;
			double E_gamma_cut = 0.3;
			float calE1 = 0;
			float calE2 = 0;
			float pcalv_e = 0;
			float pcalw_e = 0;
			byte det_gamma1 = 0;
			float pcalv_gamma1 = 0;
			float pcalw_gamma1 = 0;
			float ftcalrad_gamma1 = 0;
			byte det_gamma2 = 0;
			float pcalv_gamma2 = 0;
			float pcalw_gamma2 = 0;
			float ftcalrad_gamma2 = 0;
			LorentzVector v_e = new LorentzVector(0, 0, 0, 0);
			LorentzVector v_proton = new LorentzVector(0, 0, 0, 0);
			LorentzVector v_gamma1 = new LorentzVector(0, 0, 0, 0);
			LorentzVector v_gamma2 = new LorentzVector(0, 0, 0, 0);
			float cdchi2_proton = 100;
			short cdNDF_proton = 100;
			for(int i = 0; i < rec.rows(); i++)
			{
				int pid = rec.getInt("pid", i);
				float px = rec.getFloat("px", i);
				float py = rec.getFloat("py", i);
				float pz = rec.getFloat("pz", i);
				float vx = rec.getFloat("vx", i);
				float vy = rec.getFloat("vy", i);
				float vz = rec.getFloat("vz", i);
				byte charge = rec.getByte("charge", i);
				float beta = rec.getFloat("beta", i);
				float p = (float) Math.sqrt((px*px)+(py*py)+(pz*pz));
				if(pid == 11)
				{
					ecount++;
					if(ecount == 1)
					{
						p_e = new LorentzVector(px, py, pz, Math.sqrt((p*p)+(0.000511*0.000511)));
						v_e = new LorentzVector(vx, vy, vz, 0);
						e_pindex = i;
					}
				}
				if(pid == 2212)
				{
					protoncount++;
					if(protoncount == 1)
					{
						p_proton = new LorentzVector(px, py, pz, Math.sqrt((p*p)+(0.938272*0.938272)));
						v_proton = new LorentzVector(vx, vy, vz, 0);
						proton_pindex = i;
					}
				}
				if(pid == 22)
				{
					gammacount++;
					if(gammacount == 1)
					{	
						p_gamma1 = new LorentzVector(px, py, pz, p);
						v_gamma1 = new LorentzVector(vx, vy, vz, 0);
						gamma1_pindex = i;
					}
					if(gammacount == 2)
					{
						p_gamma2 = new LorentzVector(px, py, pz, p);
						v_gamma2 = new LorentzVector(vx, vy, vz, 0);
						gamma2_pindex = i;
					}
				}
				double SF = 0;
				byte proton_detector = 0;
				if(event.hasBank("REC::Calorimeter"))
				{
					HipoDataBank reccal = (HipoDataBank) event.getBank("REC::Calorimeter");
					for(int l = 0; l < reccal.rows(); l++)
					{
						short pindex = reccal.getShort("pindex", l);
						byte detector = reccal.getByte("detector", l);
						byte layer = reccal.getByte("layer", l);
						float lv = reccal.getFloat("lv", l);
						float lw = reccal.getFloat("lw", l);
						float energy = reccal.getFloat("energy", l);
						if(pindex == e_pindex && detector == 7 && layer == 1)
						{
							pcalv_e = lv;
							pcalw_e = lw;
						}
						if(pindex == gamma1_pindex && detector == 7 && layer == 1)
						{
							pcalv_gamma1 = lv;
							pcalw_gamma1 = lw;
							det_gamma1 = detector;
						}
						if(pindex == gamma2_pindex && detector == 7 && layer == 1)
						{
							pcalv_gamma2 = lv;
							pcalw_gamma2 = lw;
							det_gamma2 = detector;
						}
					}
				}
				if(event.hasBank("REC::ForwardTagger"))
				{
					HipoDataBank recft = (HipoDataBank) event.getBank("REC::ForwardTagger");
					for(int m = 0; m < recft.rows(); m++)
					{
						short pindex = recft.getShort("pindex", m);
						byte detector = recft.getByte("detector", m);
						float x = recft.getFloat("x", m);
						float y = recft.getFloat("y", m);
						if(pindex == gamma1_pindex && detector == 10)
						{
							ftcalrad_gamma1 = (float) Math.abs(Math.sqrt((x*x)+(y*y)));
							det_gamma1 = detector;
						}
						if(pindex == gamma2_pindex && detector == 10)
						{
							ftcalrad_gamma2 = (float) Math.abs(Math.sqrt((x*x)+(y*y)));
							det_gamma2 = detector;
						}
					}
				}
			}
			HipoDataBank rectrac = (HipoDataBank) event.getBank("REC::Track");
			byte e_detector = 0;
			byte e_sector = 0;
			byte proton_detector = 0;
			byte gamma_sector = 0;
			byte X0_sector = 0;
			for(int j = 0; j < rectrac.rows(); j++)
			{
				short pindex = rectrac.getShort("pindex", j);
				byte detector = rectrac.getByte("detector", j);
				byte sector = rectrac.getByte("sector", j);
				float chi2 = rectrac.getFloat("chi2", j);
				short NDF = rectrac.getShort("NDF", j);
				if(pindex == e_pindex)
				{
					e_detector = detector;
					e_sector = sector;
				}
				if(pindex == proton_pindex)
				{
					proton_detector = detector;
					cdchi2_proton = chi2;
					cdNDF_proton = NDF;
				}
			}
			HipoDataBank recev = (HipoDataBank) event.getBank("REC::Event");
			byte helic = recev.getByte("helicity", 0);		
			LorentzVector q = new LorentzVector(0, 0, E, Math.sqrt((E*E)+(0.000511*0.000511)));
			LorentzVector qr = new LorentzVector(0, 0, E, Math.sqrt((E*E)+(0.000511*0.000511)));
			q.sub(p_e);
			qr.sub(p_e);
			double Q2 = -1*q.mass2();
			LorentzVector X = new LorentzVector(0, 0, 0, 0.938272);
			double x_B = Q2/(2*X.e()*q.e());
			LorentzVector t = new LorentzVector(p_proton.px(), p_proton.py(), p_proton.pz(), p_proton.e());
			t.sub(X);
			X.add(q);
			q.sub(p_gamma1);
			double negt_qproton = -q.mass2();
			double W2 = X.mass2();
			double W = X.mass();
			X.sub(p_proton);
			double X_ep_m2 = X.mass2();
			LorentzVector rn_X = X;
			X.sub(p_gamma1);
			LorentzVector X_e = new LorentzVector(0, 0, E, Math.sqrt((E*E)+(0.000511*0.000511))+0.938272);
			X_e.sub(p_proton);
			X_e.sub(p_gamma1);
			LorentzVector X_proton = new LorentzVector(0, 0, E, Math.sqrt((E*E)+(0.000511*0.000511))+0.938272);
			X_proton.sub(p_e);
			LorentzVector rn_X_proton = X_proton;
			X_proton.sub(p_gamma1);
			LorentzVector X_tar = new LorentzVector(0, 0, 0, 0.938272);
			LorentzVector t_cal = new LorentzVector(X_proton.px(), X_proton.py(), X_proton.pz(), X_proton.e());
			t_cal.sub(X_tar);
			LorentzVector X_gamma1 = new LorentzVector(0, 0, E, Math.sqrt((E*E)+(0.000511*0.000511))+0.938272);
			X_gamma1.sub(p_proton);
			X_gamma1.sub(p_e);
			LorentzVector beam = new LorentzVector(0, 0, E, Math.sqrt((E*E)+(0.000511*0.000511)));
			
			double trento = ((beam.vect().cross(p_e.vect())).dot(p_proton.vect()))/Math.abs((beam.vect().cross(p_e.vect())).dot(p_proton.vect()));
			double phi = trento*57.3*Math.acos(((beam.vect().cross(p_e.vect())).dot((p_proton.vect().cross(qr.vect()))))
							/((beam.vect().cross(p_e.vect()).mag()*(p_proton.vect().cross(qr.vect())).mag())));
			if(phi < 0) phi = phi+360;
			
			int Q2xBbin = -1;
			if(x_B < 0.16) Q2xBbin = 0;
			if(x_B >= 0.16 && x_B < 0.21 && Q2 < 1.25) Q2xBbin = 1;
			if(x_B >= 0.16 && x_B < 0.21 && Q2 >= 1.25) Q2xBbin = 2;
			if(x_B >= 0.21 && x_B < 0.32 && Q2 < 1.3) Q2xBbin = 3;
			if(x_B >= 0.21 && x_B < 0.32 && Q2 >= 1.3 && Q2 < 1.65) Q2xBbin = 4;
			if(x_B >= 0.21 && x_B < 0.32 && Q2 >= 1.65) Q2xBbin = 5;
			if(x_B >= 0.32 && Q2 < 2.4) Q2xBbin = 6;
			if(x_B >= 0.32 && Q2 >= 2.4) Q2xBbin = 7;
			
			int ntbin = -1;
			if(-t_cal.mass2() < 0.2) ntbin = 0;
			if(-t_cal.mass2() >= 0.2 && -t_cal.mass2() < 0.29) ntbin = 1;
			if(-t_cal.mass2() >= 0.29 && -t_cal.mass2() < 0.38) ntbin = 2;
			if(-t_cal.mass2() >= 0.38 && -t_cal.mass2() < 0.48) ntbin = 3;
			if(-t_cal.mass2() >= 0.48 && -t_cal.mass2() < 0.6) ntbin = 4;
			if(-t_cal.mass2() >= 0.6 && -t_cal.mass2() < 0.77) ntbin = 5;
			if(-t_cal.mass2() >= 0.77 && -t_cal.mass2() < 1.04) ntbin = 6;
			if(-t_cal.mass2() >= 1.04) ntbin = 7;
			
			if(ecount == 1 && protoncount == 1 && gammacount == 1
					&& v_e.pz() > -12 && v_e.pz() < 7 && Math.abs(v_e.pz()-v_proton.pz()) < 2.5 + (2.5/p_proton.p())
					&& pcalv_e > 12 && pcalw_e > 12 && (( pcalv_gamma1 > 12 && pcalw_gamma1 > 12) || (ftcalrad_gamma1 > 8 && ftcalrad_gamma1 < 15))
					&& p_gamma1.vect().theta(p_e.vect()) > 5
					&& ((proton_detector == 6) || (proton_detector == 5 && cdNDF_proton < 10 && cdchi2_proton/cdNDF_proton < 30))
					&& 57.3*p_proton.theta() < 75)
			{
				htheta_vs_p_e.fill( p_e.p(), 57.3*p_e.theta());
				htheta_vs_phi_e_rb.fill(57.3*p_e.phi(), 57.3*p_e.theta());
				hvz_vs_p_e.fill(p_e.p(), v_e.pz());
				hvz_vs_theta_e.fill(57.3*p_e.theta(), v_e.pz());
				htheta_vs_p_proton.fill( p_proton.p(), 57.3*p_proton.theta());
				htheta_vs_phi_proton_rb.fill(57.3*p_proton.phi(), 57.3*p_proton.theta());
				hvz_vs_p_proton.fill(p_proton.p(), v_proton.pz());
				hvz_vs_theta_proton.fill(57.3*p_proton.theta(), v_proton.pz());
				htheta_vs_p_gamma.fill(p_gamma1.p(), 57.3*p_gamma1.theta());
				htheta_vs_phi_gamma_rb.fill(57.3*p_gamma1.phi(), 57.3*p_gamma1.theta());
				hvz_vs_p_gamma.fill(p_gamma1.p(), v_gamma1.pz());
				hvz_vs_theta_gamma.fill(57.3*p_gamma1.theta(), v_gamma1.pz());
				if((det_gamma1 == 7 && p_gamma1.vect().theta(X_gamma1.vect()) < 4.5 && Math.sqrt((X.px()*X.px())+(X.py()*X.py())) < 0.3
						&& X.e() < 1.25 && X_proton.mass() > 0.4 && X_proton.mass() < 1.8) && Q2 > 1 && W > 2)
				{
					htheta_cone_gamma_rb.fill(p_gamma1.vect().theta(X_gamma1.vect()));
					hX_proton_m_rb.fill(X_proton.mass());
					hX_epg_pt_rb.fill(Math.sqrt((X.px()*X.px())+(X.py()*X.py())));
					hX_epg_E_rb.fill(X.e());
					hQ2_vs_x_B_rb.fill(x_B, Q2);
					{
						
						hnegt_vs_phi_trento.fill(phi, -t_cal.mass2());
					//	int n = -1;
					//	n = (int) (phi/12);
					//	double phi_val_temp = 0;
					//	double phi_count_temp = 0;
						if(helic == 1)
						{
							hphi_CM_neghel.fill(phi);
							histGroups_phi_CM_neghel_bin.getItem(Q2xBbin).fill(phi);
							histGroups_phi_CM_neghel_tbin.getItem(Q2xBbin, ntbin ).fill(phi);
						}
						if(helic == -1)
						{	
							hphi_CM_poshel.fill(phi);
							histGroups_phi_CM_poshel_bin.getItem(Q2xBbin).fill(phi);
							histGroups_phi_CM_poshel_tbin.getItem(Q2xBbin, ntbin ).fill(phi);
						}
						/*
						if(helic == 1 || helic == -1)
						{
							phi_val_temp = phi_val.get(n);
							phi_val.set(n, phi_val_temp+phi);
							phi_count_temp = phi_count.get(n);
							phi_count.set(n, phi_count_temp+1);
						}
						*/
					}
				}
			}
			
			if(ecount == 1 && protoncount == 1 && gammacount == 2
					&& p_gamma1.vect().theta(p_gamma2.vect()) > 2
						&& v_e.pz() > -12 && v_e.pz() < 7 && Math.abs(v_e.pz()-v_proton.pz()) < 2.5 + (2.5/p_proton.p())
						&& pcalv_e > 12 && pcalw_e > 12
						&& (( pcalv_gamma1 > 12 && pcalw_gamma1 > 12) || (ftcalrad_gamma1 > 8 && ftcalrad_gamma1 < 15))
						&& (( pcalv_gamma2 > 12 && pcalw_gamma2 > 12) || (ftcalrad_gamma2 > 8 && ftcalrad_gamma2 < 15))
						&& p_gamma1.vect().theta(p_e.vect()) > 5
						&& p_gamma2.vect().theta(p_e.vect()) > 5
						&& p_gamma1.vect().theta(p_gamma2.vect()) > 2
						&& ((proton_detector == 6) || (proton_detector == 5 && cdNDF_proton < 10 && cdchi2_proton/cdNDF_proton < 30))
						&& 57.3*p_proton.theta() < 75)
			{
				LorentzVector cp_p_gamma1 = new LorentzVector (0, 0, 0, 0);
				cp_p_gamma1.copy(p_gamma1);
				X_proton.sub(p_gamma2);
				X.sub(p_gamma2);
				p_gamma1.add(p_gamma2);
				if(((det_gamma1 == 7 && det_gamma2 == 7
						&& p_gamma1.vect().theta(X_gamma1.vect()) < 4.5
						&& Math.sqrt((X.px()*X.px())+(X.py()*X.py())) < 0.3
						&& X.e() < 1.25 && X_proton.mass() > 0.4 && X_proton.mass() < 1.8) && Q2 > 1 && W > 2) && Q2 > 1 && W > 2)
				{
					
					hX_proton_m_pi0_rb.fill(X_proton.mass());
					hX_epgg_pt_rb.fill(Math.sqrt((X.px()*X.px())+(X.py()*X.py())));
					hX_epgg_E_rb.fill(X.e());
					htheta_cone_pi0_rb.fill(p_gamma1.vect().theta(X_gamma1.vect()));
					if(p_gamma1.mass() > 0.101 && p_gamma1.mass() < 0.167)
					{
						h2gamma_m.fill(p_gamma1.mass());
						hX_2gamma_m.fill(X_gamma1.mass());
						
						if(helic == 1)
						{
							hphi_CM_neghel_pi0.fill(phi);
							histGroups_phi_CM_neghel_pi0_bin.getItem(Q2xBbin).fill(phi);
							histGroups_phi_CM_neghel_pi0_tbin.getItem(Q2xBbin, ntbin).fill(phi);
						}
						if(helic == -1)
						{
							hphi_CM_poshel_pi0.fill(phi);
							histGroups_phi_CM_poshel_pi0_bin.getItem(Q2xBbin).fill(phi);
							histGroups_phi_CM_poshel_pi0_tbin.getItem(Q2xBbin, ntbin).fill(phi);
						}
					}
				}
			}			
		}
	}
	
	static void processEvent_pi0_sim(DataEvent event, int eventCounter, List<Double> phi_val, List<Double> phi_count)
	{
		if(event.hasBank("REC::Particle") && event.hasBank("REC::Track"))
		{
			HipoDataBank rec = (HipoDataBank) event.getBank("REC::Particle");
			int ecount = 0;
			int protoncount = 0;
			int gammacount = 0;
			int e_pindex = -1;
			int proton_pindex = -1;
			int gamma1_pindex = -1;
			int gamma2_pindex = -1;
			LorentzVector p_e = new LorentzVector(0, 0, 0, 0.000511);
			LorentzVector p_proton = new LorentzVector(0, 0, 0, 0.938272);
			LorentzVector p_gamma1 = new LorentzVector(0, 0, 0, 0);
			LorentzVector p_gamma2 = new LorentzVector(0, 0, 0, 0);
			LorentzVector p_X0 = new LorentzVector(0, 0, 0, 0);
			double E = 7.546;
			double E_gamma_cut = 0.3;
			float calE1 = 0;
			float calE2 = 0;
			float pcalv_e = 0;
			float pcalw_e = 0;
			byte det_gamma1 = 0;
			float pcalv_gamma1 = 0;
			float pcalw_gamma1 = 0;
			float ftcalrad_gamma1 = 0;
			byte det_gamma2 = 0;
			float pcalv_gamma2 = 0;
			float pcalw_gamma2 = 0;
			float ftcalrad_gamma2 = 0;
			byte calsector = 0;
			LorentzVector v_e = new LorentzVector(0, 0, 0, 0);
			LorentzVector v_proton = new LorentzVector(0, 0, 0, 0);
			LorentzVector v_gamma1 = new LorentzVector(0, 0, 0, 0);
			LorentzVector v_gamma2 = new LorentzVector(0, 0, 0, 0);
			float cdchi2_proton = 100;
			short cdNDF_proton = 100;
			for(int i = 0; i < rec.rows(); i++)
			{
				int pid = rec.getInt("pid", i);
				float px = rec.getFloat("px", i);
				float py = rec.getFloat("py", i);
				float pz = rec.getFloat("pz", i);
				float vx = rec.getFloat("vx", i);
				float vy = rec.getFloat("vy", i);
				float vz = rec.getFloat("vz", i);
				byte charge = rec.getByte("charge", i);
				float beta = rec.getFloat("beta", i);
				float p = (float) Math.sqrt((px*px)+(py*py)+(pz*pz));
				if(pid == 11)
				{
					ecount++;
					if(ecount == 1)
					{
						p_e = new LorentzVector(px, py, pz, Math.sqrt((p*p)+(0.000511*0.000511)));
						v_e = new LorentzVector(vx, vy, vz, 0);
						e_pindex = i;
					}
				}
				if(pid == 2212)
				{
					protoncount++;
					if(protoncount == 1)
					{
						p_proton = new LorentzVector(px, py, pz, Math.sqrt((p*p)+(0.938272*0.938272)));
						v_proton = new LorentzVector(vx, vy, vz, 0);
						proton_pindex = i;
					}
				}
				if(pid == 22)
				{
					gammacount++;
					if(gammacount == 1)
					{	
					//	p_gamma1 = new LorentzVector(px+(0.25*px/p), py+(0.25*py/p), pz+(0.25*pz/p), p+0.25);
						p_gamma1 = new LorentzVector(px, py, pz, p);
						v_gamma1 = new LorentzVector(vx, vy, vz, 0);
						gamma1_pindex = i;
					}
					if(gammacount == 2)
					{
					//	p_gamma2 = new LorentzVector(px+(0.25*px/p), py+(0.25*py/p), pz+(0.25*pz/p), p+0.25);
						p_gamma2 = new LorentzVector(px, py, pz, p);
						v_gamma2 = new LorentzVector(vx, vy, vz, 0);
						gamma2_pindex = i;
					}
				}
				byte proton_detector = 0;
				if(event.hasBank("REC::Calorimeter"))
				{
					HipoDataBank reccal = (HipoDataBank) event.getBank("REC::Calorimeter");
					for(int l = 0; l < reccal.rows(); l++)
					{
						short pindex = reccal.getShort("pindex", l);
						byte detector = reccal.getByte("detector", l);
						byte layer = reccal.getByte("layer", l);
						float lv = reccal.getFloat("lv", l);
						float lw = reccal.getFloat("lw", l);
						float x = reccal.getFloat("x", l);
						float y = reccal.getFloat("y", l);
						if(pindex == e_pindex && detector == 7 && layer == 1)
						{
							pcalv_e = lv;
							pcalw_e = lw;
						}
						if(pindex == gamma1_pindex && detector == 7 && layer == 1)
						{
							pcalv_gamma1 = lv;
							pcalw_gamma1 = lw;
							det_gamma1 = detector;
						}
						if(pindex == gamma2_pindex && detector == 7 && layer == 1)
						{
							pcalv_gamma2 = lv;
							pcalw_gamma2 = lw;
							det_gamma2 = detector;
						}
					}
				}
				if(event.hasBank("REC::ForwardTagger"))
				{
					HipoDataBank recft = (HipoDataBank) event.getBank("REC::ForwardTagger");
					for(int m = 0; m < recft.rows(); m++)
					{
						short pindex = recft.getShort("pindex", m);
						byte detector = recft.getByte("detector", m);
						float x = recft.getFloat("x", m);
						float y = recft.getFloat("y", m);
						if(pindex == gamma1_pindex && detector == 10)
						{
							ftcalrad_gamma1 = (float) Math.abs(Math.sqrt((x*x)+(y*y)));
							det_gamma1 = detector;
						}
						if(pindex == gamma2_pindex && detector == 10)
						{
							ftcalrad_gamma2 = (float) Math.abs(Math.sqrt((x*x)+(y*y)));
							det_gamma2 = detector;
						}
					}
				}
			}
			HipoDataBank rectrac = (HipoDataBank) event.getBank("REC::Track");
			byte e_detector = 0;
			byte e_sector = 0;
			byte proton_detector = 0;
			byte gamma_sector = 0;
			byte X0_sector = 0;
			for(int j = 0; j < rectrac.rows(); j++)
			{
				short pindex = rectrac.getShort("pindex", j);
				byte detector = rectrac.getByte("detector", j);
				byte sector = rectrac.getByte("sector", j);
				float chi2 = rectrac.getFloat("chi2", j);
				short NDF = rectrac.getShort("NDF", j);
				if(pindex == e_pindex)
				{
					e_detector = detector;
					e_sector = sector;
				}
				if(pindex == proton_pindex)
				{
					proton_detector = detector;
					cdchi2_proton = chi2;
					cdNDF_proton = NDF;
				}
			}
			HipoDataBank recev = (HipoDataBank) event.getBank("REC::Event");
			byte helic = recev.getByte("helicity", 0);		
			LorentzVector q = new LorentzVector(0, 0, E, Math.sqrt((E*E)+(0.000511*0.000511)));
			LorentzVector qr = new LorentzVector(0, 0, E, Math.sqrt((E*E)+(0.000511*0.000511)));
			q.sub(p_e);
			qr.sub(p_e);
			double Q2 = -1*q.mass2();
			LorentzVector X = new LorentzVector(0, 0, 0, 0.938272);
			double x_B = Q2/(2*X.e()*q.e());
			LorentzVector t = new LorentzVector(p_proton.px(), p_proton.py(), p_proton.pz(), p_proton.e());
			t.sub(X);
			X.add(q);
			q.sub(p_gamma1);
			double negt_qproton = -q.mass2();
			double W2 = X.mass2();
			double W = X.mass();
			X.sub(p_proton);
			double X_ep_m2 = X.mass2();
			LorentzVector rn_X1 = X;
			LorentzVector rn_X2 = X;
			X.sub(p_gamma1);
			LorentzVector X_e = new LorentzVector(0, 0, E, Math.sqrt((E*E)+(0.000511*0.000511))+0.938272);
			X_e.sub(p_proton);
			X_e.sub(p_gamma1);
			LorentzVector X_proton = new LorentzVector(0, 0, E, Math.sqrt((E*E)+(0.000511*0.000511))+0.938272);
			X_proton.sub(p_e);
			LorentzVector rn_X_proton1 = X_proton;
			LorentzVector rn_X_proton2 = X_proton;
			X_proton.sub(p_gamma1);
			LorentzVector X_tar = new LorentzVector(0, 0, 0, 0.938272);
			LorentzVector t_cal = new LorentzVector(X_proton.px(), X_proton.py(), X_proton.pz(), X_proton.e());
			t_cal.sub(X_tar);
			LorentzVector X_gamma1 = new LorentzVector(0, 0, E, Math.sqrt((E*E)+(0.000511*0.000511))+0.938272);
			X_gamma1.sub(p_proton);
			X_gamma1.sub(p_e);
			LorentzVector beam = new LorentzVector(0, 0, E, Math.sqrt((E*E)+(0.000511*0.000511)));
			
			double trento = ((beam.vect().cross(p_e.vect())).dot(p_proton.vect()))/Math.abs((beam.vect().cross(p_e.vect())).dot(p_proton.vect()));
			double phi = trento*57.3*Math.acos(((beam.vect().cross(p_e.vect())).dot((p_proton.vect().cross(qr.vect()))))
							/((beam.vect().cross(p_e.vect()).mag()*(p_proton.vect().cross(qr.vect())).mag())));
			if(phi < 0) phi = phi+360;
			
			int Q2xBbin = -1;
			if(x_B < 0.16) Q2xBbin = 0;
			if(x_B >= 0.16 && x_B < 0.21 && Q2 < 1.25) Q2xBbin = 1;
			if(x_B >= 0.16 && x_B < 0.21 && Q2 >= 1.25) Q2xBbin = 2;
			if(x_B >= 0.21 && x_B < 0.32 && Q2 < 1.3) Q2xBbin = 3;
			if(x_B >= 0.21 && x_B < 0.32 && Q2 >= 1.3 && Q2 < 1.65) Q2xBbin = 4;
			if(x_B >= 0.21 && x_B < 0.32 && Q2 >= 1.65) Q2xBbin = 5;
			if(x_B >= 0.32 && Q2 < 2.4) Q2xBbin = 6;
			if(x_B >= 0.32 && Q2 >= 2.4) Q2xBbin = 7;
			
			int ntbin = -1;
			if(-t_cal.mass2() < 0.2) ntbin = 0;
			if(-t_cal.mass2() >= 0.2 && -t_cal.mass2() < 0.29) ntbin = 1;
			if(-t_cal.mass2() >= 0.29 && -t_cal.mass2() < 0.38) ntbin = 2;
			if(-t_cal.mass2() >= 0.38 && -t_cal.mass2() < 0.48) ntbin = 3;
			if(-t_cal.mass2() >= 0.48 && -t_cal.mass2() < 0.6) ntbin = 4;
			if(-t_cal.mass2() >= 0.6 && -t_cal.mass2() < 0.77) ntbin = 5;
			if(-t_cal.mass2() >= 0.77 && -t_cal.mass2() < 1.04) ntbin = 6;
			if(-t_cal.mass2() >= 1.04) ntbin = 7;
			
			if(ecount == 1 && protoncount == 1 && gammacount == 1
					&& v_e.pz() > -12 && v_e.pz() < 7 && Math.abs(v_e.pz()-v_proton.pz()) < 2.5 + (2.5/p_proton.p())
					&& pcalv_e > 12 && pcalw_e > 12 && (( pcalv_gamma1 > 12 && pcalw_gamma1 > 12)
					|| (57.3*p_gamma1.theta() < 4.5 && 57.3*p_gamma1.theta() > 2.5))
					&& p_gamma1.vect().theta(p_e.vect()) > 5
					&& ((proton_detector == 6) || (proton_detector == 5 && cdNDF_proton < 10 && cdchi2_proton/cdNDF_proton < 30))
					&& 57.3*p_proton.theta() < 75)
			{
				/*
				if(pcalv_e > 12) System.out.println("pcalv_e = " + pcalv_e);
				if(ftcalrad_gamma1 > 8 && ftcalrad_gamma1 < 15) System.out.println("ftcalrad_gamma1 = " + ftcalrad_gamma1);
				if(cdNDF_proton < 10) System.out.println("cdNDF_proton = " + cdNDF_proton);
				if(cdchi2_proton/cdNDF_proton < 30) System.out.println("cdchi2_proton/cdNDF_proton = " + cdchi2_proton/cdNDF_proton);
				*/
				
				if((det_gamma1 == 7 && p_gamma1.vect().theta(X_gamma1.vect()) < 4.5 && Math.sqrt((X.px()*X.px())+(X.py()*X.py())) < 0.3
						&& X.e() < 1.25 && X_proton.mass() > 0.4 && X_proton.mass() < 1.8) && Q2 > 1 && W > 2)
				{
					{
					//	int n = -1;
					//	n = (int) (phi/12);
					//	double phi_val_temp = 0;
					//	double phi_count_temp = 0;
					//	if(helic == 1 || helic == -1)
						{
							
							hphi_CM_1gamma_pi0_sim.fill(phi);
							histGroups_phi_CM_1gamma_pi0_sim_bin.getItem(Q2xBbin).fill(phi);
							histGroups_phi_CM_1gamma_pi0_sim_tbin.getItem(Q2xBbin, ntbin).fill(phi);
						}
						/*
						if(helic == 1 || helic == -1)
						{
							phi_val_temp = phi_val.get(n);
							phi_val.set(n, phi_val_temp+phi);
							phi_count_temp = phi_count.get(n);
							phi_count.set(n, phi_count_temp+1);
						}
						*/
					}
				}
			}
			
			if(ecount == 1 && protoncount == 1 && gammacount == 2
						&& v_e.pz() > -12 && v_e.pz() < 7 && Math.abs(v_e.pz()-v_proton.pz()) < 2.5 + (2.5/p_proton.p())
						&& pcalv_e > 12 && pcalw_e > 12
						&& (( pcalv_gamma1 > 12 && pcalw_gamma1 > 12) || (57.3*p_gamma1.theta() < 4.5 && 57.3*p_gamma1.theta() > 2.5))
						&& (( pcalv_gamma2 > 12 && pcalw_gamma2 > 12) || (57.3*p_gamma2.theta() < 4.5 && 57.3*p_gamma2.theta() > 2.5))
						&& p_gamma1.vect().theta(p_e.vect()) > 5
						&& p_gamma2.vect().theta(p_e.vect()) > 5
						&& p_gamma1.vect().theta(p_gamma2.vect()) > 2
						&& ((proton_detector == 6) || (proton_detector == 5 && cdNDF_proton < 10 && cdchi2_proton/cdNDF_proton < 30))
						&& 57.3*p_proton.theta() < 75)
			{
				LorentzVector cp_p_gamma1 = new LorentzVector (0, 0, 0, 0);
				cp_p_gamma1.copy(p_gamma1);
				X_proton.sub(p_gamma2);
				X.sub(p_gamma2);
				p_gamma1.add(p_gamma2);
				if(((det_gamma1 == 7 && det_gamma2 == 7
						&& p_gamma1.vect().theta(X_gamma1.vect()) < 4.5
						&& Math.sqrt((X.px()*X.px())+(X.py()*X.py())) < 0.3
						&& X.e() < 1.25 && X_proton.mass() > 0.4 && X_proton.mass() < 1.8) && Q2 > 1 && W > 2) && Q2 > 1 && W > 2)
				{
					if(p_gamma1.mass() > 0.101 && p_gamma1.mass() < 0.167)
					{	
						double x1 = 612.5*cp_p_gamma1.px()/cp_p_gamma1.pz();
						double y1 = 612.5*cp_p_gamma1.py()/cp_p_gamma1.pz();
						double phi1 = 57.3*cp_p_gamma1.phi();
						double theta1 = 57.3*cp_p_gamma1.theta();
						double x2 = 612.5*p_gamma2.px()/p_gamma2.pz();
						double y2 = 612.5*p_gamma2.py()/p_gamma2.pz();
						double phi2 = 57.3*p_gamma2.phi();
						double theta2 = 57.3*p_gamma2.theta();
					//	if(helic == 1 || helic == -1)
						{
							hphi_CM_2gamma_pi0_sim.fill(phi);
							histGroups_phi_CM_2gamma_pi0_sim_bin.getItem(Q2xBbin).fill(phi);
							histGroups_phi_CM_2gamma_pi0_sim_tbin.getItem(Q2xBbin, ntbin).fill(phi);
						}
					}
				}
			}
		}
	}
	
	public static void main(String[] args) {
		
		phi_CM_poshel_bin_histos();
		phi_CM_neghel_bin_histos();
		phi_CM_poshel_tbin_histos();
		phi_CM_neghel_tbin_histos();
		phi_CM_poshel_pi0_bin_histos();
		phi_CM_neghel_pi0_bin_histos();
		phi_CM_poshel_pi0_tbin_histos();
		phi_CM_neghel_pi0_tbin_histos();
		phi_CM_1gamma_pi0_sim_bin_histos();
		phi_CM_2gamma_pi0_sim_bin_histos();
		phi_CM_1gamma_pi0_sim_tbin_histos();
		phi_CM_2gamma_pi0_sim_tbin_histos();
		
		reader.open("C:/Users/joshtanj/Documents/download/merged_skim16_6535MeV.hipo");
		reader_1gamma_pi0_sim.open("C:/Users/joshtanj/Documents/download/merged_simu_rgk_exc_1photon_dvpi0p_7.5GeV.hipo");
		reader_pi0_sim.open("C:/Users/joshtanj/Documents/download/merged_simu_rgk_exc_dvpi0p_7.5GeV.hipo");
		
		List<Double> phi_val = new ArrayList<>();
		List<Double> phi_count = new ArrayList<>();
				
		for(int pent = 0; pent < 30; pent++)
		{
			phi_val.add((double) 0);
			phi_count.add((double) 0);
		}
		
		int eventCounter = 1;
		while(reader.hasEvent())// && eventCounter < 10000000)
		{
			processEvent(reader.getNextEvent(), eventCounter, phi_val, phi_count);
			if(eventCounter%50000 == 0) System.out.println("Event: " + eventCounter);
			eventCounter++;
		}
		int Nevent = eventCounter-1;
		System.out.println("Number of events: " + Nevent);
		
		int eventCounter_1gamma_sim = 1;
		while(reader_1gamma_pi0_sim.hasEvent())// && eventCounter_1gamma_sim < 1000)
		{
			processEvent_pi0_sim(reader_1gamma_pi0_sim.getNextEvent(), eventCounter_1gamma_sim, phi_val, phi_count);
			if(eventCounter_1gamma_sim%50000 == 0) System.out.println("#pi^0(1#gamma)_sim event: " + eventCounter_1gamma_sim);
			eventCounter_1gamma_sim++;
		}
		int Nevent_1gamma = eventCounter_1gamma_sim-1;
		System.out.println("Number of #pi^0(1#gamma)_sim events: " + Nevent_1gamma);
		
		int eventCounter_pi0_sim = 1;
		while(reader_pi0_sim.hasEvent())// && eventCounter_pi0_sim < 1000)
		{
			processEvent_pi0_sim(reader_pi0_sim.getNextEvent(), eventCounter_pi0_sim, phi_val, phi_count);
			if(eventCounter_pi0_sim%50000 == 0) System.out.println("#pi^0_sim event: " + eventCounter_pi0_sim);
			eventCounter_pi0_sim++;
		}
		int Nevent_pi0 = eventCounter_pi0_sim-1;
		System.out.println("Number of #pi^0_sim events: " + Nevent_pi0);
		
		List<Double> bsa = new ArrayList<>();
		List<Double> pi0_bsa = new ArrayList<>();
		List<Double> pi0_sim_r = new ArrayList<>();
		List<Double> pi0_R = new ArrayList<>();
		List<Double> pi0_sim_f = new ArrayList<>();
		List<Double> corr_pi0_sim_bsa = new ArrayList<>();
		
		for(int i = 0; i < 30; i++)
		{
			bsa.add((1/0.863)*(hphi_CM_poshel.getBinContent(i)-hphi_CM_neghel.getBinContent(i))
						/(hphi_CM_poshel.getBinContent(i)+hphi_CM_neghel.getBinContent(i)));
			pi0_bsa.add((1/0.863)*(hphi_CM_poshel_pi0.getBinContent(i)-hphi_CM_neghel_pi0.getBinContent(i))
							/(hphi_CM_poshel_pi0.getBinContent(i)+hphi_CM_neghel_pi0.getBinContent(i)));
			pi0_sim_r.add(hphi_CM_1gamma_pi0_sim.getBinContent(i)/hphi_CM_2gamma_pi0_sim.getBinContent(i));
			pi0_R.add((hphi_CM_poshel_pi0.getBinContent(i)+hphi_CM_neghel_pi0.getBinContent(i))
							/((hphi_CM_poshel.getBinContent(i)+hphi_CM_neghel.getBinContent(i))));
			pi0_sim_f.add((hphi_CM_poshel_pi0.getBinContent(i)+hphi_CM_neghel_pi0.getBinContent(i))*pi0_sim_r.get(i)
							/(hphi_CM_poshel.getBinContent(i)+hphi_CM_neghel.getBinContent(i)));
			corr_pi0_sim_bsa.add((bsa.get(i)-(pi0_sim_f.get(i)*pi0_bsa.get(i)))/(1-pi0_sim_f.get(i)));
		}
		IndexedList<Double> bsa_bin = new IndexedList<>(2);
		IndexedList<Double> pi0_bsa_bin = new IndexedList<>(2);
		IndexedList<Double> pi0_sim_r_bin = new IndexedList<>(2);
		IndexedList<Double> pi0_R_bin = new IndexedList<>(2);
		IndexedList<Double> pi0_sim_f_bin = new IndexedList<>(2);
		IndexedList<Double> corr_pi0_sim_bsa_bin = new IndexedList<>(2);
		for(int ibin = 0; ibin  < 8; ibin++)
		{
			for(int ient = 0; ient < 30; ient++)
			{
				bsa_bin.add((1/0.863)*(histGroups_phi_CM_poshel_bin.getItem(ibin).getBinContent(ient)
								-histGroups_phi_CM_neghel_bin.getItem(ibin).getBinContent(ient))
								/(histGroups_phi_CM_poshel_bin.getItem(ibin).getBinContent(ient)
								+histGroups_phi_CM_neghel_bin.getItem(ibin).getBinContent(ient)), ibin, ient);
				pi0_bsa_bin.add((1/0.863)*(histGroups_phi_CM_poshel_pi0_bin.getItem(ibin).getBinContent(ient)
									-histGroups_phi_CM_neghel_pi0_bin.getItem(ibin).getBinContent(ient))
									/(histGroups_phi_CM_poshel_pi0_bin.getItem(ibin).getBinContent(ient)
									+histGroups_phi_CM_neghel_pi0_bin.getItem(ibin).getBinContent(ient)), ibin, ient);
				pi0_sim_r_bin.add(histGroups_phi_CM_1gamma_pi0_sim_bin.getItem(ibin).getBinContent(ient)
								/histGroups_phi_CM_2gamma_pi0_sim_bin.getItem(ibin).getBinContent(ient), ibin, ient);
				pi0_R_bin.add((histGroups_phi_CM_poshel_pi0_bin.getItem(ibin).getBinContent(ient)
								+histGroups_phi_CM_neghel_pi0_bin.getItem(ibin).getBinContent(ient))
								/(histGroups_phi_CM_poshel_bin.getItem(ibin).getBinContent(ient)
								+histGroups_phi_CM_neghel_bin.getItem(ibin).getBinContent(ient)), ibin, ient);
				pi0_sim_f_bin.add((histGroups_phi_CM_poshel_pi0_bin.getItem(ibin).getBinContent(ient)
									+histGroups_phi_CM_neghel_pi0_bin.getItem(ibin).getBinContent(ient))*pi0_sim_r_bin.getItem(ibin, ient)
									/(histGroups_phi_CM_poshel_bin.getItem(ibin).getBinContent(ient)
									+histGroups_phi_CM_neghel_bin.getItem(ibin).getBinContent(ient)), ibin, ient);
				corr_pi0_sim_bsa_bin.add((bsa_bin.getItem(ibin, ient)-(pi0_sim_f_bin.getItem(ibin, ient)*pi0_bsa_bin.getItem(ibin, ient)))
											/(1-pi0_sim_f_bin.getItem(ibin, ient)), ibin, ient);
			}
		}
		List<Double> bsaerr = new ArrayList<>();
		IndexedList<Double> bsaerr_bin = new IndexedList<>(2);
		List<Double> pi0_bsaerr = new ArrayList<>();
		IndexedList<Double> pi0_bsaerr_bin = new IndexedList<>(2);
		List<Double> pi0_sim_rerr = new ArrayList<>();
		IndexedList<Double> pi0_sim_rerr_bin = new IndexedList<>(2);
		List<Double> pi0_sim_ferr = new ArrayList<>();
		IndexedList<Double> pi0_sim_ferr_bin = new IndexedList<>(2);
		List<Double> corr_pi0_sim_bsaerr = new ArrayList<>();
		IndexedList<Double> corr_pi0_sim_bsaerr_bin = new IndexedList<>(2);
		
		for(int j = 0; j < 30; j++)
		{
			bsaerr.add((1/0.863)*Math.sqrt((1-(0.863*0.863*bsa.get(j)*bsa.get(j)))/(hphi_CM_poshel.getBinContent(j)+hphi_CM_neghel.getBinContent(j))));
			pi0_bsaerr.add((1/0.863)*Math.sqrt((1-(0.863*0.863*pi0_bsa.get(j)*pi0_bsa.get(j)))
							/(hphi_CM_poshel_pi0.getBinContent(j)+hphi_CM_neghel_pi0.getBinContent(j))));
			pi0_sim_rerr.add(Math.sqrt(pi0_sim_r.get(j)*(1+pi0_sim_r.get(j))/hphi_CM_2gamma_pi0_sim.getBinContent(j)));
			pi0_sim_ferr.add(pi0_sim_r.get(j)*Math.sqrt(pi0_R.get(j)*(1+pi0_R.get(j))
							/(hphi_CM_poshel.getBinContent(j)+hphi_CM_neghel.getBinContent(j)))+pi0_R.get(j)*pi0_sim_rerr.get(j));
			corr_pi0_sim_bsaerr.add((1/0.863)*2*(1/(hphi_CM_poshel.getBinContent(j)+hphi_CM_neghel.getBinContent(j)))
									*(1/(hphi_CM_poshel.getBinContent(j)+hphi_CM_neghel.getBinContent(j)))
									*Math.sqrt((hphi_CM_neghel.getBinContent(j)*hphi_CM_neghel.getBinContent(j)
									*(hphi_CM_poshel.getBinContent(j)+hphi_CM_poshel_pi0.getBinContent(j)*hphi_CM_poshel_pi0.getBinContent(j)
									*pi0_sim_rerr.get(j)*pi0_sim_rerr.get(j)+pi0_sim_r.get(j)*pi0_sim_r.get(j)*hphi_CM_poshel_pi0.getBinContent(j)))
									+(hphi_CM_poshel.getBinContent(j)*hphi_CM_poshel.getBinContent(j)
									*(hphi_CM_neghel.getBinContent(j)+hphi_CM_neghel_pi0.getBinContent(j)*hphi_CM_neghel_pi0.getBinContent(j)
									*pi0_sim_rerr.get(j)*pi0_sim_rerr.get(j)+pi0_sim_r.get(j)*pi0_sim_r.get(j)*hphi_CM_neghel_pi0.getBinContent(j)))));
		}
		for(int jbin = 0; jbin  < 8; jbin++)
		{
			for(int jent = 0; jent < 30; jent++)
			{
				bsaerr_bin.add((1/0.863)*Math.sqrt((1-(0.863*0.863*bsa_bin.getItem(jbin, jent)*bsa_bin.getItem(jbin, jent)))
								/(histGroups_phi_CM_poshel_bin.getItem(jbin).getBinContent(jent)
								+histGroups_phi_CM_neghel_bin.getItem(jbin).getBinContent(jent))), jbin, jent);
				pi0_bsaerr_bin.add((1/0.863)*Math.sqrt((1-(0.863*0.863*pi0_bsa_bin.getItem(jbin, jent)*pi0_bsa_bin.getItem(jbin, jent)))
										/(histGroups_phi_CM_poshel_pi0_bin.getItem(jbin).getBinContent(jent)
										+histGroups_phi_CM_neghel_pi0_bin.getItem(jbin).getBinContent(jent))), jbin, jent);
				pi0_sim_rerr_bin.add(Math.sqrt(pi0_sim_r_bin.getItem(jbin, jent)*(1+pi0_sim_r_bin.getItem(jbin, jent))
										/histGroups_phi_CM_2gamma_pi0_sim_bin.getItem(jbin).getBinContent(jent)), jbin, jent);
				pi0_sim_ferr_bin.add(pi0_sim_r_bin.getItem(jbin, jent)*Math.sqrt(pi0_R_bin.getItem(jbin, jent)*(1+pi0_R_bin.getItem(jbin, jent))
										/(histGroups_phi_CM_poshel_bin.getItem(jbin).getBinContent(jent)
										+histGroups_phi_CM_neghel_bin.getItem(jbin).getBinContent(jent)))
										+pi0_R_bin.getItem(jbin, jent)*pi0_sim_rerr_bin.getItem(jbin, jent), jbin, jent);
				corr_pi0_sim_bsaerr_bin.add((1/0.863)*2*(1/(histGroups_phi_CM_poshel_bin.getItem(jbin).getBinContent(jent)
												+histGroups_phi_CM_neghel_bin.getItem(jbin).getBinContent(jent)))
												*(1/(histGroups_phi_CM_poshel_bin.getItem(jbin).getBinContent(jent)
												+histGroups_phi_CM_neghel_bin.getItem(jbin).getBinContent(jent)))
												*Math.sqrt((histGroups_phi_CM_neghel_bin.getItem(jbin).getBinContent(jent)
												*histGroups_phi_CM_neghel_bin.getItem(jbin).getBinContent(jent)
												*(histGroups_phi_CM_poshel_bin.getItem(jbin).getBinContent(jent)
												+histGroups_phi_CM_poshel_pi0_bin.getItem(jbin).getBinContent(jent)
												*histGroups_phi_CM_poshel_pi0_bin.getItem(jbin).getBinContent(jent)
												*pi0_sim_rerr_bin.getItem(jbin, jent)*pi0_sim_rerr_bin.getItem(jbin, jent)
												+pi0_sim_r_bin.getItem(jbin, jent)*pi0_sim_r_bin.getItem(jbin, jent)
												*histGroups_phi_CM_poshel_pi0_bin.getItem(jbin).getBinContent(jent)))
												+(histGroups_phi_CM_poshel_bin.getItem(jbin).getBinContent(jent)
												*histGroups_phi_CM_poshel_bin.getItem(jbin).getBinContent(jent)
												*(histGroups_phi_CM_neghel_bin.getItem(jbin).getBinContent(jent)
												+histGroups_phi_CM_neghel_pi0_bin.getItem(jbin).getBinContent(jent)
												*histGroups_phi_CM_neghel_pi0_bin.getItem(jbin).getBinContent(jent)
												*pi0_sim_rerr_bin.getItem(jbin, jent)*pi0_sim_rerr_bin.getItem(jbin, jent)
												+pi0_sim_r_bin.getItem(jbin, jent)*pi0_sim_r_bin.getItem(jbin, jent)
												*histGroups_phi_CM_neghel_pi0_bin.getItem(jbin).getBinContent(jent)))), jbin, jent);
			}
		}
		
		IndexedList<Double> bsa_tbin = new IndexedList<>(3);
		IndexedList<Double> pi0_bsa_tbin = new IndexedList<>(3);
		IndexedList<Double> pi0_sim_r_tbin = new IndexedList<>(3);
		IndexedList<Double> pi0_R_tbin = new IndexedList<>(3);
		IndexedList<Double> pi0_sim_f_tbin = new IndexedList<>(3);
		IndexedList<Double> corr_pi0_sim_bsa_tbin = new IndexedList<>(3);
		for(int ibin = 0; ibin  < 8; ibin++)
		{
			for (int itbin = 0; itbin < 8; itbin++)
			{
				for(int ient = 0; ient < 30; ient++)
				{
					bsa_tbin.add((1/0.863)*(histGroups_phi_CM_poshel_tbin.getItem(ibin, itbin).getBinContent(ient)
									-histGroups_phi_CM_neghel_tbin.getItem(ibin, itbin).getBinContent(ient))
									/((histGroups_phi_CM_poshel_tbin.getItem(ibin, itbin).getBinContent(ient)
									+histGroups_phi_CM_neghel_tbin.getItem(ibin, itbin).getBinContent(ient))), ibin, itbin, ient);
					pi0_bsa_tbin.add((1/0.863)*(histGroups_phi_CM_poshel_pi0_tbin.getItem(ibin, itbin).getBinContent(ient)
										-histGroups_phi_CM_neghel_pi0_tbin.getItem(ibin, itbin).getBinContent(ient))
										/(histGroups_phi_CM_poshel_pi0_tbin.getItem(ibin, itbin).getBinContent(ient)
										+histGroups_phi_CM_neghel_pi0_tbin.getItem(ibin, itbin).getBinContent(ient)), ibin, itbin, ient);
					pi0_sim_r_tbin.add(histGroups_phi_CM_1gamma_pi0_sim_tbin.getItem(ibin, itbin).getBinContent(ient)
										/histGroups_phi_CM_2gamma_pi0_sim_tbin.getItem(ibin, itbin).getBinContent(ient), ibin, itbin, ient);
					pi0_R_tbin.add((histGroups_phi_CM_poshel_pi0_tbin.getItem(ibin, itbin).getBinContent(ient)
									+histGroups_phi_CM_neghel_pi0_tbin.getItem(ibin, itbin).getBinContent(ient))
									/(histGroups_phi_CM_poshel_tbin.getItem(ibin, itbin).getBinContent(ient)
									+histGroups_phi_CM_neghel_tbin.getItem(ibin, itbin).getBinContent(ient)), ibin, itbin, ient);
					pi0_sim_f_tbin.add((histGroups_phi_CM_poshel_pi0_tbin.getItem(ibin, itbin).getBinContent(ient)
											+histGroups_phi_CM_neghel_pi0_tbin.getItem(ibin, itbin).getBinContent(ient))*pi0_sim_r_tbin.getItem(ibin, itbin, ient)
											/(histGroups_phi_CM_poshel_tbin.getItem(ibin, itbin).getBinContent(ient)
											+histGroups_phi_CM_neghel_tbin.getItem(ibin, itbin).getBinContent(ient)), ibin, itbin, ient);
					corr_pi0_sim_bsa_tbin.add((bsa_tbin.getItem(ibin, itbin, ient)-(pi0_sim_f_tbin.getItem(ibin, itbin, ient)*pi0_bsa_tbin.getItem(ibin, itbin, ient)))
												/(1-pi0_sim_f_tbin.getItem(ibin, itbin, ient)), ibin, itbin, ient);
				}
			}
		}
		IndexedList<Double> bsaerr_tbin = new IndexedList<>(3);
		IndexedList<Double> pi0_bsaerr_tbin = new IndexedList<>(3);
		IndexedList<Double> pi0_sim_rerr_tbin = new IndexedList<>(3);
		IndexedList<Double> pi0_sim_ferr_tbin = new IndexedList<>(3);
		IndexedList<Double> corr_pi0_sim_bsaerr_tbin = new IndexedList<>(3);
		for(int jbin = 0; jbin  < 8; jbin++)
		{
			for (int jtbin = 0; jtbin < 8; jtbin++)
			{
				for(int jent = 0; jent < 30; jent++)
				{
					bsaerr_tbin.add((1/0.863)*Math.sqrt((1-(0.863*0.863*bsa_tbin.getItem(jbin, jtbin, jent)*bsa_tbin.getItem(jbin, jtbin, jent)))
									/(histGroups_phi_CM_poshel_tbin.getItem(jbin, jtbin).getBinContent(jent)
									+histGroups_phi_CM_neghel_tbin.getItem(jbin, jtbin).getBinContent(jent))), jbin, jtbin, jent);
					pi0_bsaerr_tbin.add((1/0.863)*Math.sqrt((1-(0.863*0.863*pi0_bsa_tbin.getItem(jbin, jtbin, jent)*pi0_bsa_tbin.getItem(jbin, jtbin, jent)))
											/(histGroups_phi_CM_poshel_pi0_tbin.getItem(jbin, jtbin).getBinContent(jent)
											+histGroups_phi_CM_neghel_pi0_tbin.getItem(jbin, jtbin).getBinContent(jent))), jbin, jtbin, jent);
					pi0_sim_rerr_tbin.add(Math.sqrt(pi0_sim_r_tbin.getItem(jbin, jtbin, jent)*(1+pi0_sim_r_tbin.getItem(jbin, jtbin, jent))
											/histGroups_phi_CM_2gamma_pi0_sim_tbin.getItem(jbin, jtbin).getBinContent(jent)), jbin, jtbin, jent);
					pi0_sim_ferr_tbin.add(pi0_sim_r_tbin.getItem(jbin, jtbin, jent)
											*Math.sqrt(pi0_R_tbin.getItem(jbin, jtbin, jent)*(1+pi0_R_tbin.getItem(jbin, jtbin, jent))
											/(histGroups_phi_CM_poshel_tbin.getItem(jbin, jtbin).getBinContent(jent)
											+histGroups_phi_CM_neghel_tbin.getItem(jbin, jtbin).getBinContent(jent)))
											+pi0_R_tbin.getItem(jbin, jtbin, jent)*pi0_sim_rerr_tbin.getItem(jbin, jtbin, jent), jbin, jtbin, jent);
					corr_pi0_sim_bsaerr_tbin.add((1/0.863)*2*(1/(histGroups_phi_CM_poshel_tbin.getItem(jbin, jtbin).getBinContent(jent)
													+histGroups_phi_CM_neghel_tbin.getItem(jbin, jtbin).getBinContent(jent)))
													*(1/(histGroups_phi_CM_poshel_tbin.getItem(jbin, jtbin).getBinContent(jent)
													+histGroups_phi_CM_neghel_tbin.getItem(jbin, jtbin).getBinContent(jent)))
													*Math.sqrt((histGroups_phi_CM_neghel_tbin.getItem(jbin, jtbin).getBinContent(jent)
													*histGroups_phi_CM_neghel_tbin.getItem(jbin, jtbin).getBinContent(jent)
													*(histGroups_phi_CM_poshel_tbin.getItem(jbin, jtbin).getBinContent(jent)
													+histGroups_phi_CM_poshel_pi0_tbin.getItem(jbin, jtbin).getBinContent(jent)
													*histGroups_phi_CM_poshel_pi0_tbin.getItem(jbin, jtbin).getBinContent(jent)
													*pi0_sim_rerr_tbin.getItem(jbin, jtbin, jent)*pi0_sim_rerr_tbin.getItem(jbin, jtbin, jent)
													+pi0_sim_r_tbin.getItem(jbin, jtbin, jent)*pi0_sim_r_tbin.getItem(jbin, jtbin, jent)
													*histGroups_phi_CM_poshel_pi0_tbin.getItem(jbin, jtbin).getBinContent(jent)))
													+(histGroups_phi_CM_poshel_tbin.getItem(jbin, jtbin).getBinContent(jent)
													*histGroups_phi_CM_poshel_tbin.getItem(jbin, jtbin).getBinContent(jent)
													*(histGroups_phi_CM_neghel_tbin.getItem(jbin, jtbin).getBinContent(jent)
													+histGroups_phi_CM_neghel_pi0_tbin.getItem(jbin, jtbin).getBinContent(jent)
													*histGroups_phi_CM_neghel_pi0_tbin.getItem(jbin, jtbin).getBinContent(jent)
													*pi0_sim_rerr_tbin.getItem(jbin, jtbin, jent)*pi0_sim_rerr_tbin.getItem(jbin, jtbin, jent)
													+pi0_sim_r_tbin.getItem(jbin, jtbin, jent)*pi0_sim_r_tbin.getItem(jbin, jtbin, jent)
													*histGroups_phi_CM_neghel_pi0_tbin.getItem(jbin, jtbin).getBinContent(jent)))), jbin, jtbin, jent);
				}
			}
		}
		
		double[] phibc = new double[] {6, 18, 30, 42, 54, 66, 78, 90, 102, 114, 126, 138, 150, 162, 174,
				 						186, 198, 210, 222, 234, 246, 258, 270, 282, 294, 306, 318, 330, 342, 354};
		GraphErrors hBSAy = new GraphErrors();
		GraphErrors pi0_hBSAy = new GraphErrors();
		GraphErrors pi0_sim_hr = new GraphErrors();
		GraphErrors pi0_sim_hf = new GraphErrors();
		GraphErrors corr_pi0_sim_hBSAy = new GraphErrors();
		
	//	System.out.println("==========================================================================================");
	//	System.out.println("ep#gamma BSA: x #deltax y #deltay");
		for(int k = 0; k < 30; k++)
		{	
	//		System.out.println(phibc[k] + "	" + 0 + "	" + bsa.get(k) + "	" + bsaerr.get(k));
			hBSAy.addPoint(phibc[k], bsa.get(k), 0, bsaerr.get(k));
		}
	//	System.out.println("==========================================================================================");
	//	System.out.println("#pi^0 BSA: x #deltax y #deltay");
		for(int k = 0; k < 30; k++)
		{	
	//		System.out.println(phibc[k] + "	" + 0 + "	" + pi0_bsa.get(k) + "	" + pi0_bsaerr.get(k));
			pi0_hBSAy.addPoint(phibc[k], pi0_bsa.get(k), 0, pi0_bsaerr.get(k));
		}
	//	System.out.println("==========================================================================================");
	//	System.out.println("r: x #deltax y #deltay");
		for(int k = 0; k < 30; k++)
		{	
	//		System.out.println(phibc[k] + "	" + 0 + "	" + pi0_sim_r.get(k) + "	" + pi0_sim_rerr.get(k));
			pi0_sim_hr.addPoint(phibc[k], pi0_sim_r.get(k), 0, pi0_sim_rerr.get(k));
		}
	//	System.out.println("==========================================================================================");
	//	System.out.println("f: x #deltax y #deltay");
		for(int k = 0; k < 30; k++)
		{	
	//		System.out.println(phibc[k] + "	" + 0 + "	" + pi0_sim_f.get(k) + "	" + pi0_sim_ferr.get(k));
			pi0_sim_hf.addPoint(phibc[k], pi0_sim_f.get(k), 0, pi0_sim_ferr.get(k));
		}
	//	System.out.println("==========================================================================================");
	//	System.out.println("#pi^0 Corrected BSA: x #deltax y #deltay");
		for(int k = 0; k < 30; k++)
		{	
	//		System.out.println(phibc[k] + "	" + 0 + "	" + corr_pi0_sim_bsa.get(k) + "	" + corr_pi0_sim_bsaerr.get(k));
			corr_pi0_sim_hBSAy.addPoint(phibc[k], corr_pi0_sim_bsa.get(k), 0, corr_pi0_sim_bsaerr.get(k));
		}
		
		IndexedList<GraphErrors> hBSAy_bin = new IndexedList<GraphErrors>(1);
		IndexedList<GraphErrors> pi0_hBSAy_bin = new IndexedList<GraphErrors>(1);
		IndexedList<GraphErrors> pi0_sim_hr_bin = new IndexedList<GraphErrors>(1);
		IndexedList<GraphErrors> pi0_sim_hf_bin = new IndexedList<GraphErrors>(1);
		IndexedList<GraphErrors> corr_pi0_sim_hBSAy_bin = new IndexedList<GraphErrors>(1);
	//	System.out.println("==========================================================================================");
	//	System.out.println("BSA: x #deltax y #deltay");
		for(int kbin = 0; kbin  < 8; kbin++)
		{
	//		System.out.println("Bin No.: " + kbin);
			GraphErrors hBSAy_bin_ini = new GraphErrors();
			hBSAy_bin.add(hBSAy_bin_ini, kbin);
			for(int kent = 0; kent < 30; kent++)
			{
				if(bsa_bin.hasItem(kbin, kent) && bsaerr_bin.hasItem(kbin, kent))
				{
	//				System.out.println(phibc[kent] + "	" + 0 + "	" + bsa_bin.getItem(kbin, kent) + "	" + bsaerr_bin.getItem(kbin, kent));
					hBSAy_bin.getItem(kbin).addPoint(phibc[kent], bsa_bin.getItem(kbin, kent), 0, bsaerr_bin.getItem(kbin, kent));
				}
			}
		}
	//	System.out.println("==========================================================================================");
	//	System.out.println("#pi^0 BSA: x #deltax y #deltay");
		for(int kbin = 0; kbin  < 8; kbin++)
		{
	//		System.out.println("Bin No.: " + kbin);
			GraphErrors pi0_hBSAy_bin_ini = new GraphErrors();
			pi0_hBSAy_bin.add(pi0_hBSAy_bin_ini, kbin);
			for(int kent = 0; kent < 30; kent++)
			{
				if(pi0_bsa_bin.hasItem(kbin, kent) && pi0_bsaerr_bin.hasItem(kbin, kent))
				{
	//				System.out.println(phibc[kent] + "	" + 0 + "	" + pi0_bsa_bin.getItem(kbin, kent) + "	" + pi0_bsaerr_bin.getItem(kbin, kent));
					pi0_hBSAy_bin.getItem(kbin).addPoint(phibc[kent], pi0_bsa_bin.getItem(kbin, kent), 0, pi0_bsaerr_bin.getItem(kbin, kent));
				}
			}
		}
	//	System.out.println("==========================================================================================");
	//	System.out.println("r: x #deltax y #deltay");
		for(int kbin = 0; kbin  < 8; kbin++)
		{
	//		System.out.println("Bin No.: " + kbin);
			GraphErrors pi0_sim_hr_bin_ini = new GraphErrors();
			pi0_sim_hr_bin.add(pi0_sim_hr_bin_ini, kbin);
			for(int kent = 0; kent < 30; kent++)
			{
				if(pi0_sim_r_bin.hasItem(kbin, kent) && pi0_sim_rerr_bin.hasItem(kbin, kent))
				{
	//				System.out.println(phibc[kent] + "	" + 0 + "	" + pi0_sim_r_bin.getItem(kbin, kent) + "	" + pi0_sim_rerr_bin.getItem(kbin, kent));
					pi0_sim_hr_bin.getItem(kbin).addPoint(phibc[kent], pi0_sim_r_bin.getItem(kbin, kent), 0, pi0_sim_rerr_bin.getItem(kbin, kent));
				}
			}
		}
	//	System.out.println("==========================================================================================");
	//	System.out.println("f: x #deltax y #deltay");
		for(int kbin = 0; kbin  < 8; kbin++)
		{
	//		System.out.println("Bin No.: " + kbin);
			GraphErrors pi0_sim_hf_bin_ini = new GraphErrors();
			pi0_sim_hf_bin.add(pi0_sim_hf_bin_ini, kbin);
			for(int kent = 0; kent < 30; kent++)
			{
				if(pi0_sim_f_bin.hasItem(kbin, kent))
				{
	//				System.out.println(phibc[kent] + "	" + 0 + "	" + pi0_sim_f_bin.getItem(kbin, kent) + "	" + pi0_sim_ferr_bin.getItem(kbin, kent));
					pi0_sim_hf_bin.getItem(kbin).addPoint(phibc[kent], pi0_sim_f_bin.getItem(kbin, kent), 0, pi0_sim_ferr_bin.getItem(kbin, kent));
				}
			}
		}
	//	System.out.println("==========================================================================================");
	//	System.out.println("#pi^0 Corrected BSA: x #deltax y #deltay");
		for(int kbin = 0; kbin  < 8; kbin++)
		{
	//		System.out.println("Bin No.: " + kbin);
			GraphErrors corr_pi0_sim_hBSAy_bin_ini = new GraphErrors();
			corr_pi0_sim_hBSAy_bin.add(corr_pi0_sim_hBSAy_bin_ini, kbin);
			for(int kent = 0; kent < 30; kent++)
			{
				if(corr_pi0_sim_bsa_bin.hasItem(kbin, kent) && corr_pi0_sim_bsaerr_bin.hasItem(kbin, kent))
				{
	//				System.out.println(phibc[kent] + "	" + 0 + "	" + corr_pi0_sim_bsa_bin.getItem(kbin, kent) + "	"
	//									+ corr_pi0_sim_bsaerr_bin.getItem(kbin, kent));
					corr_pi0_sim_hBSAy_bin.getItem(kbin).addPoint(phibc[kent], corr_pi0_sim_bsa_bin.getItem(kbin, kent),
																	0, corr_pi0_sim_bsaerr_bin.getItem(kbin, kent));
				}
			}
		}
		IndexedList<GraphErrors> hBSAy_tbin = new IndexedList<GraphErrors>(2);
		IndexedList<GraphErrors> halpha_vs_t = new IndexedList<GraphErrors>(1);
		IndexedList<Double> beta_vs_t = new IndexedList<>(2);
		IndexedList<Double> betaerr_vs_t = new IndexedList<>(2);
		IndexedList<GraphErrors> pi0_hBSAy_tbin = new IndexedList<GraphErrors>(2);
		IndexedList<GraphErrors> pi0_halpha_vs_t = new IndexedList<GraphErrors>(1);
		IndexedList<GraphErrors> pi0_sim_hf_tbin = new IndexedList<GraphErrors>(2);
		IndexedList<GraphErrors> corr_pi0_sim_hBSAy_tbin = new IndexedList<GraphErrors>(2);
		IndexedList<GraphErrors> corr_pi0_sim_halpha_vs_t = new IndexedList<GraphErrors>(1);
		double[] tbc = new double[] {0.14, 0.25, 0.34, 0.44, 0.54, 0.68, 0.91, 1.47};
	//	System.out.println("==========================================================================================");
	//	System.out.println("BSA: #alpha vs -t");
		for(int kbin = 0; kbin  < 8; kbin++)
		{
	//		System.out.println("	Q^2-x_B Bin No.: " + kbin);
			GraphErrors halpha_vs_t_ini = new GraphErrors();
			halpha_vs_t.add(halpha_vs_t_ini, kbin);
			for(int ktbin = 0; ktbin < 8; ktbin++)
			{
	//			System.out.println("		-t Bin No.: " + ktbin);
				GraphErrors hBSAy_tbin_ini = new GraphErrors();
				hBSAy_tbin.add(hBSAy_tbin_ini, kbin, ktbin);
				for(int kent = 0; kent < 30; kent++)
				{
					if(bsa_tbin.hasItem(kbin, ktbin, kent) && bsaerr_tbin.hasItem(kbin, ktbin, kent) && bsaerr_tbin.getItem(kbin, ktbin, kent) < 0.5)
					{
	//					System.out.println("			" + phibc[kent] + "	" + 0 + "	"
	//							+ bsa_tbin.getItem(kbin, ktbin, kent) + "	" + bsaerr_tbin.getItem(kbin, ktbin, kent));
						hBSAy_tbin.getItem(kbin, ktbin).addPoint(phibc[kent], bsa_tbin.getItem(kbin, ktbin, kent), 0, bsaerr_tbin.getItem(kbin, ktbin, kent));
					}
				}
				F1D fBSA= new F1D("fBSA", "[A]*sin(x/57.3)/(1+[B]*cos(x/57.3))", 0, 360);
				fBSA.setParameter(0, 0.1);
	//			fBSA.setParLimits(0, 0, 1.0);
				fBSA.setParameter(1, -0.5);
				DataFitter.fit(fBSA, hBSAy_tbin.getItem(kbin, ktbin), "Q");
				if(fBSA.parameter(0).value() > 0.099 && fBSA.parameter(0).value() <0.101)
				{
					fBSA.setParameter(0, 0.25);
	//				fBSA.setParLimits(0, 0, 1.0);
					fBSA.setParameter(1, -0.5);
					DataFitter.fit(fBSA, hBSAy_tbin.getItem(kbin, ktbin), "Q");
				}
				if(fBSA.parameter(0).value() > 0.24 && fBSA.parameter(0).value() < 0.26)
				{
					fBSA.setParameter(0, 0.0);
	//				fBSA.setParLimits(0, 0, 1.0);
					fBSA.setParameter(1, -0.5);
					DataFitter.fit(fBSA, hBSAy_tbin.getItem(kbin, ktbin), "Q");
				}
				halpha_vs_t.getItem(kbin).addPoint(tbc[ktbin], fBSA.parameter(0).value(), 0, fBSA.parameter(0).error());
				beta_vs_t.add(fBSA.parameter(1).value(), kbin, ktbin);
				betaerr_vs_t.add(fBSA.parameter(1).error(),  kbin, ktbin);
			}
		}
	//	System.out.println("==========================================================================================");
	//	System.out.println("#pi^0 BSA: #alpha vs -t");
		for(int kbin = 0; kbin  < 8; kbin++)
		{
	//		System.out.println("	Q^2-x_B Bin No.: " + kbin);
			GraphErrors pi0_halpha_vs_t_ini = new GraphErrors();
			pi0_halpha_vs_t.add(pi0_halpha_vs_t_ini, kbin);
			for(int ktbin = 0; ktbin < 8; ktbin++)
			{
	//			System.out.println("		-t Bin No.: " + ktbin);
				GraphErrors pi0_hBSAy_tbin_ini = new GraphErrors();
				pi0_hBSAy_tbin.add(pi0_hBSAy_tbin_ini, kbin, ktbin);
				for(int kent = 0; kent < 30; kent++)
				{
					if(pi0_bsa_tbin.hasItem(kbin, ktbin, kent) && pi0_bsaerr_tbin.hasItem(kbin, ktbin, kent) && pi0_bsaerr_tbin.getItem(kbin, ktbin, kent) < 0.5)
					{
	//					System.out.println("			" + phibc[kent] + "	" + 0 + "	"
	//							+ pi0_bsa_tbin.getItem(kbin, ktbin, kent) + "	" + pi0_bsaerr_tbin.getItem(kbin, ktbin, kent));
						pi0_hBSAy_tbin.getItem(kbin, ktbin).addPoint(phibc[kent], pi0_bsa_tbin.getItem(kbin, ktbin, kent),
																		0, pi0_bsaerr_tbin.getItem(kbin, ktbin, kent));
					}
				}
				F1D fBSA= new F1D("fBSA", "[A]*sin(x/57.3)", 0, 360);
				fBSA.setParameter(0, 0.2);
				fBSA.setParLimits(0, 0, 0.5);
				DataFitter.fit(fBSA, pi0_hBSAy_tbin.getItem(kbin, ktbin), "Q");
				pi0_halpha_vs_t.getItem(kbin).addPoint(tbc[ktbin], fBSA.parameter(0).value(), 0, fBSA.parameter(0).error());
			}
		}
		
	//	System.out.println("==========================================================================================");
	//	System.out.println("#pi^0 Corrected BSA: #alpha vs -t");
		for(int kbin = 0; kbin  < 8; kbin++)
		{
	//		System.out.println("	Q^2-x_B Bin No.: " + kbin);
			GraphErrors corr_pi0_sim_halpha_vs_t_ini = new GraphErrors();
			corr_pi0_sim_halpha_vs_t.add(corr_pi0_sim_halpha_vs_t_ini, kbin);
			for(int ktbin = 0; ktbin < 8; ktbin++)
			{
	//			System.out.println("-t Bin No.: " + ktbin);
				GraphErrors corr_pi0_sim_hBSAy_tbin_ini = new GraphErrors();
				corr_pi0_sim_hBSAy_tbin.add(corr_pi0_sim_hBSAy_tbin_ini, kbin, ktbin);
				for(int kent = 0; kent < 30; kent++)
				{
					if(corr_pi0_sim_bsa_tbin.hasItem(kbin, ktbin, kent) && corr_pi0_sim_bsaerr_tbin.hasItem(kbin, ktbin, kent)
							&& corr_pi0_sim_bsaerr_tbin.getItem(kbin, ktbin, kent) < 0.5)
					{
	//					System.out.println("			" + phibc[kent] + "	" + 0 + "	"
	//							+ corr_pi0_sim_bsa_tbin.getItem(kbin, ktbin, kent) + "	" + corr_pi0_sim_bsaerr_tbin.getItem(kbin, ktbin, kent));
						corr_pi0_sim_hBSAy_tbin.getItem(kbin, ktbin).addPoint(phibc[kent], corr_pi0_sim_bsa_tbin.getItem(kbin, ktbin, kent),
																		0, corr_pi0_sim_bsaerr_tbin.getItem(kbin, ktbin, kent));
					}
				}
				F1D fBSA= new F1D("fBSA", "[A]*sin(x/57.3)/(1+[B]*cos(x/57.3))", 0, 360);
				fBSA.setParameter(0, 0.2);
				fBSA.setParameter(1, beta_vs_t.getItem(kbin, ktbin));
				fBSA.setParLimits(1, beta_vs_t.getItem(kbin, ktbin)+0.000001,
									beta_vs_t.getItem(kbin, ktbin)-0.000001);
				DataFitter.fit(fBSA, corr_pi0_sim_hBSAy_tbin.getItem(kbin, ktbin), "Q");
				/*
				if(fBSA.parameter(0).value() > 0.099 && fBSA.parameter(0).value() <0.101)
				{
					fBSA.setParameter(0, 0.25);
					fBSA.setParLimits(0, -0.05, 0.5);
					fBSA.setParameter(1, beta_vs_t.getItem(kbin, ktbin));
					fBSA.setParLimits(1, beta_vs_t.getItem(kbin, ktbin)-betaerr_vs_t.getItem(kbin, ktbin),
										beta_vs_t.getItem(kbin, ktbin)+betaerr_vs_t.getItem(kbin, ktbin));
					DataFitter.fit(fBSA, corr_pi0_sim_hBSAy_tbin.getItem(kbin, ktbin), "Q");
				}
				if(fBSA.parameter(0).value() > 0.24 && fBSA.parameter(0).value() < 0.26)
				{
					fBSA.setParameter(0, 0.0);
					fBSA.setParLimits(0, -0.05, 0.5);
					fBSA.setParameter(1, beta_vs_t.getItem(kbin, ktbin));
					fBSA.setParLimits(1, beta_vs_t.getItem(kbin, ktbin)-betaerr_vs_t.getItem(kbin, ktbin),
										beta_vs_t.getItem(kbin, ktbin)+betaerr_vs_t.getItem(kbin, ktbin));
					DataFitter.fit(fBSA, corr_pi0_sim_hBSAy_tbin.getItem(kbin, ktbin), "Q");
				}
				*/
	//			System.out.println("	" + tbc[ktbin] + "	" + 0 + "	" + fBSA.parameter(0).value() + "	" + fBSA.parameter(0).error());
				corr_pi0_sim_halpha_vs_t.getItem(kbin).addPoint(tbc[ktbin], fBSA.parameter(0).value(), 0, fBSA.parameter(0).error());
			}
		}
		/*
		System.out.println("==========================================================================================");
		System.out.println("#pi^0 Corrected BSA: #alpha vs -t:");
		for(int kbin = 0; kbin  < 8; kbin++)
		{
			System.out.println("	Q^2-x_B Bin No.: " + kbin);
			GraphErrors corr_pi0_sim_halpha_vs_t_ini = new GraphErrors();
			corr_pi0_sim_halpha_vs_t.add(corr_pi0_sim_halpha_vs_t_ini, kbin);
			for(int ktbin = 0; ktbin < 8; ktbin++)
			{
				System.out.println("		-t Bin No.: " + ktbin);
				GraphErrors corr_pi0_sim_hBSAy_tbin_ini = new GraphErrors();
				corr_pi0_sim_hBSAy_tbin.add(corr_pi0_sim_hBSAy_tbin_ini, kbin, ktbin);
				GraphErrors pi0_sim_hf_tbin_ini = new GraphErrors();
				pi0_sim_hf_tbin.add(pi0_sim_hf_tbin_ini, kbin, ktbin);
				double f = 0;
				double ferr = 0;
				double corr_pi0_halpha = 0;
				double corr_pi0_halphaerr = 0;
				for(int kent = 0; kent < 30; kent++)
				{
					if(pi0_sim_f_tbin.hasItem(kbin, ktbin, kent) && pi0_sim_ferr_tbin.hasItem(kbin, ktbin, kent))
					{
	//					System.out.println("			" + phibc[kent] + "	" + 0 + "	"
	//							+ pi0_sim_f_tbin.getItem(kbin, ktbin, kent) + "	" + pi0_sim_ferr_tbin.getItem(kbin, ktbin, kent));
						pi0_sim_hf_tbin.getItem(kbin, ktbin).addPoint(phibc[kent], pi0_sim_f_tbin.getItem(kbin, ktbin, kent),
																		0, pi0_sim_ferr_tbin.getItem(kbin, ktbin, kent));
					}
				}
				if(halpha_vs_t.getItem(kbin).getDataX(ktbin) > 0 && pi0_halpha_vs_t.getItem(kbin).getDataX(ktbin) > 0)
				{
					F1D ff_sim= new F1D("ff_sim", "[a]*(exp(-0.5((x-[m1])/[s])((x-[m1])/[s]))+exp(-0.5((x-[m2])/[s])((x-[m2])/[s]))) ", 0, 360);
					ff_sim.setParameter(0, 0.5);
					ff_sim.setParLimits(0, 0, 1);
					ff_sim.setParameter(1, 90);
					ff_sim.setParLimits(1, 65, 115);
					ff_sim.setParameter(2, 45);
					ff_sim.setParLimits(2, 10, 80);
					ff_sim.setParameter(3, 270);
					ff_sim.setParLimits(3, 245, 395);
					DataFitter.fit(ff_sim, pi0_sim_hf_tbin.getItem(kbin, ktbin), "Q");
					double a = ff_sim.parameter(0).value();
					double da = ff_sim.parameter(0).error();
					double m1 = ff_sim.parameter(1).value();
					double dm1 = ff_sim.parameter(1).error();
					double s = ff_sim.parameter(2).value();
					double ds = ff_sim.parameter(2).error();
					double m2 = ff_sim.parameter(3).value();
					double dm2 = ff_sim.parameter(3).error();
					double e1 = Math.exp(-0.5*((90-m1)/s)*((90-m1)/s));
					double e2 = Math.exp(-0.5*((90-m2)/s)*((90-m2)/s));
					f = a*(Math.exp(-0.5*((90-m1)/s)*((90-m1)/s))+Math.exp(-0.5*((90-m2)/s)*((90-m2)/s)));
					ferr = Math.sqrt(((e1+e2)*da)*((e1+e2)*da)+(a*((90-m1)/s)*e1*dm1)*(a*((90-m1)/s)*e1*dm1)+(a*((90-m2)/s)*e2*dm2)*(a*((90-m2)/s)*e2*dm2)
										+(a*(((90-m1)/s)*e1+((90-m2)/s)*e2)*ds/s)*(a*(((90-m1)/s)*e1+((90-m2)/s)*e2)*ds/s));
					corr_pi0_halpha = ((halpha_vs_t.getItem(kbin).getDataY(ktbin)-f*pi0_halpha_vs_t.getItem(kbin).getDataY(ktbin))/(1-f));
					corr_pi0_halphaerr = Math.sqrt((halpha_vs_t.getItem(kbin).getDataEY(ktbin)/(1-f))
													*(halpha_vs_t.getItem(kbin).getDataEY(ktbin)/(1-f))
													+(pi0_halpha_vs_t.getItem(kbin).getDataEY(ktbin)/(1-f))
													*(pi0_halpha_vs_t.getItem(kbin).getDataEY(ktbin)/(1-f))
													+((halpha_vs_t.getItem(kbin).getDataY(ktbin)-pi0_halpha_vs_t.getItem(kbin).getDataY(ktbin))*ferr/((1-f)*(1-f)))
													*((halpha_vs_t.getItem(kbin).getDataY(ktbin)-pi0_halpha_vs_t.getItem(kbin).getDataY(ktbin))*ferr/((1-f)*(1-f))));
					corr_pi0_sim_halpha_vs_t.getItem(kbin).addPoint(tbc[ktbin], corr_pi0_halpha, 0, corr_pi0_halphaerr);
					System.out.println(halpha_vs_t.getItem(kbin).getDataY(ktbin) + "	" + halpha_vs_t.getItem(kbin).getDataEY(ktbin) + "	"
											+ pi0_halpha_vs_t.getItem(kbin).getDataY(ktbin) + "	" + pi0_halpha_vs_t.getItem(kbin).getDataEY(ktbin)+ "	"
											+ f + "	" + ferr);
				}
			}
		}
		*/
		JFrame framepi0 = new JFrame("#gamma#gamma Mass");
		framepi0.setSize(1000, 500);
		EmbeddedCanvas canpi0 = new EmbeddedCanvas();
		framepi0.add(canpi0);
		framepi0.setLocationRelativeTo(null);
		framepi0.setVisible(true);
		canpi0.divide(2, 1);
		canpi0.cd(0);
		canpi0.setFont("Arial");
		h2gamma_m.setTitle("#gamma#gamma Mass");
		h2gamma_m.setTitleX("Mass [GeV]");
		h2gamma_m.setTitleY("Counts");
		h2gamma_m.setOptStat(10);
		canpi0.getPad(0).setTitleFontSize(32);
		canpi0.getPad(0).setAxisTitleFontSize(32);
		canpi0.getPad(0).setAxisLabelFontSize(24);
		canpi0.getPad(0).setStatBoxFontSize(18);
		canpi0.draw(h2gamma_m, "same");
		canpi0.cd(1);
		canpi0.setFont("Arial");
		hX_2gamma_m.setTitle("sqrt((p_X_(ep->e'p'))^2)");
		hX_2gamma_m.setTitleX("sqrt((p_X)^2) [GeV]");
		hX_2gamma_m.setTitleY("Counts");
		hX_2gamma_m.setOptStat(10);
		canpi0.getPad(1).setTitleFontSize(32);
		canpi0.getPad(1).setAxisTitleFontSize(32);
		canpi0.getPad(1).setAxisLabelFontSize(24);
		canpi0.getPad(1).setStatBoxFontSize(18);
		canpi0.draw(hX_2gamma_m, "same");
		
		JFrame framee = new JFrame("e^- Kinematics");
		framee.setSize(1000, 1000);
		EmbeddedCanvas cane = new EmbeddedCanvas();
		framee.add(cane);
		framee.setLocationRelativeTo(null);
		framee.setVisible(true);
		cane.divide(2, 2);
		cane.cd(0);
		cane.setFont("Arial");
		htheta_vs_p_e.setTitle("e^- #theta vs. p");
		htheta_vs_p_e.setTitleX("p [GeV]");
		htheta_vs_p_e.setTitleY("#theta [deg]");
		cane.getPad(0).setTitleFontSize(32);
		cane.getPad(0).setAxisTitleFontSize(32);
		cane.getPad(0).setAxisLabelFontSize(24);
		cane.getPad(0).setStatBoxFontSize(18);
		cane.draw(htheta_vs_p_e, "same");
		cane.cd(1);
		cane.setFont("Arial");
		htheta_vs_phi_e_rb.setTitle("e^- #theta vs. #phi");
		htheta_vs_phi_e_rb.setTitleX("#phi [deg]");
		htheta_vs_phi_e_rb.setTitleY("#theta [deg]");
		cane.getPad(1).setTitleFontSize(32);
		cane.getPad(1).setAxisTitleFontSize(32);
		cane.getPad(1).setAxisLabelFontSize(24);
		cane.getPad(1).setStatBoxFontSize(18);
		cane.draw(htheta_vs_phi_e_rb, "same");
		cane.cd(2);
		cane.setFont("Arial");
		hvz_vs_p_e.setTitle("e^- v_z vs. p");
		hvz_vs_p_e.setTitleX("p [GeV]");
		hvz_vs_p_e.setTitleY("v_z [cm]");
		cane.getPad(2).setTitleFontSize(32);
		cane.getPad(2).setAxisTitleFontSize(32);
		cane.getPad(2).setAxisLabelFontSize(24);
		cane.getPad(2).setStatBoxFontSize(18);
		cane.draw(hvz_vs_p_e, "same");
		cane.cd(3);
		cane.setFont("Arial");
		hvz_vs_theta_e.setTitle("e^- v_z vs. #theta");
		hvz_vs_theta_e.setTitleX("#theta [deg]");
		hvz_vs_theta_e.setTitleY("v_z [cm]");
		cane.getPad(3).setTitleFontSize(32);
		cane.getPad(3).setAxisTitleFontSize(32);
		cane.getPad(3).setAxisLabelFontSize(24);
		cane.getPad(3).setStatBoxFontSize(18);
		cane.draw(hvz_vs_theta_e, "same");
		
		JFrame frameproton = new JFrame("proton Kinematics");
		frameproton.setSize(1000, 1000);
		EmbeddedCanvas canproton = new EmbeddedCanvas();
		frameproton.add(canproton);
		frameproton.setLocationRelativeTo(null);
		frameproton.setVisible(true);
		canproton.divide(2, 2);
		canproton.cd(0);
		canproton.setFont("Arial");
		htheta_vs_p_proton.setTitle("proton #theta vs. p");
		htheta_vs_p_proton.setTitleX("p [GeV]");
		htheta_vs_p_proton.setTitleY("#theta [deg]");
		canproton.getPad(0).setTitleFontSize(32);
		canproton.getPad(0).setAxisTitleFontSize(32);
		canproton.getPad(0).setAxisLabelFontSize(24);
		canproton.getPad(0).setStatBoxFontSize(18);
		canproton.draw(htheta_vs_p_proton, "same");
		canproton.cd(1);
		canproton.setFont("Arial");
		htheta_vs_phi_proton_rb.setTitle("proton #theta vs. #phi");
		htheta_vs_phi_proton_rb.setTitleX("#phi [deg]");
		htheta_vs_phi_proton_rb.setTitleY("#theta [deg]");
		canproton.getPad(1).setTitleFontSize(32);
		canproton.getPad(1).setAxisTitleFontSize(32);
		canproton.getPad(1).setAxisLabelFontSize(24);
		canproton.getPad(1).setStatBoxFontSize(18);
		canproton.draw(htheta_vs_phi_proton_rb, "same");
		canproton.cd(2);
		canproton.setFont("Arial");
		hvz_vs_p_proton.setTitle("proton v_z vs. p");
		hvz_vs_p_proton.setTitleX("p [GeV]");
		hvz_vs_p_proton.setTitleY("v_z [cm]");
		canproton.getPad(2).setTitleFontSize(32);
		canproton.getPad(2).setAxisTitleFontSize(32);
		canproton.getPad(2).setAxisLabelFontSize(24);
		canproton.getPad(2).setStatBoxFontSize(18);
		canproton.draw(hvz_vs_p_proton, "same");
		canproton.cd(3);
		canproton.setFont("Arial");
		hvz_vs_theta_proton.setTitle("proton v_z vs. #theta");
		hvz_vs_theta_proton.setTitleX("#theta [deg]");
		hvz_vs_theta_proton.setTitleY("v_z [cm]");
		canproton.getPad(3).setTitleFontSize(32);
		canproton.getPad(3).setAxisTitleFontSize(32);
		canproton.getPad(3).setAxisLabelFontSize(24);
		canproton.getPad(3).setStatBoxFontSize(18);
		canproton.draw(hvz_vs_theta_proton, "same");
		
		JFrame framegamma = new JFrame("#gamma Kinematics");
		framegamma.setSize(1000, 1000);
		EmbeddedCanvas cangamma = new EmbeddedCanvas();
		framegamma.add(cangamma);
		framegamma.setLocationRelativeTo(null);
		framegamma.setVisible(true);
		cangamma.divide(2, 2);
		cangamma.cd(0);
		cangamma.setFont("Arial");
		htheta_vs_p_gamma.setTitle("#gamma #theta vs. p");
		htheta_vs_p_gamma.setTitleX("p [GeV]");
		htheta_vs_p_gamma.setTitleY("#theta [deg]");
		cangamma.getPad(0).setTitleFontSize(32);
		cangamma.getPad(0).setAxisTitleFontSize(32);
		cangamma.getPad(0).setAxisLabelFontSize(24);
		cangamma.getPad(0).setStatBoxFontSize(18);
		cangamma.draw(htheta_vs_p_gamma, "same");
		cangamma.cd(1);
		cangamma.setFont("Arial");
		htheta_vs_phi_gamma_rb.setTitle("#gamma #theta vs. #phi");
		htheta_vs_phi_gamma_rb.setTitleX("#phi [deg]");
		htheta_vs_phi_gamma_rb.setTitleY("#theta [deg]");
		cangamma.getPad(1).setTitleFontSize(32);
		cangamma.getPad(1).setAxisTitleFontSize(32);
		cangamma.getPad(1).setAxisLabelFontSize(24);
		cangamma.getPad(1).setStatBoxFontSize(18);
		cangamma.draw(htheta_vs_phi_gamma_rb, "same");
		cangamma.cd(2);
		cangamma.setFont("Arial");
		hvz_vs_p_gamma.setTitle("#gamma v_z vs. p");
		hvz_vs_p_gamma.setTitleX("p [GeV]");
		hvz_vs_p_gamma.setTitleY("v_z [cm]");
		cangamma.getPad(2).setTitleFontSize(32);
		cangamma.getPad(2).setAxisTitleFontSize(32);
		cangamma.getPad(2).setAxisLabelFontSize(24);
		cangamma.getPad(2).setStatBoxFontSize(18);
		cangamma.draw(hvz_vs_p_gamma, "same");
		cangamma.cd(3);
		cangamma.setFont("Arial");
		hvz_vs_theta_gamma.setTitle("#gamma v_z vs. #theta");
		hvz_vs_theta_gamma.setTitleX("#theta [deg]");
		hvz_vs_theta_gamma.setTitleY("v_z [cm]");
		cangamma.getPad(3).setTitleFontSize(32);
		cangamma.getPad(3).setAxisTitleFontSize(32);
		cangamma.getPad(3).setAxisLabelFontSize(24);
		cangamma.getPad(3).setStatBoxFontSize(18);
		cangamma.draw(hvz_vs_theta_gamma, "same");
		
		JFrame frameX = new JFrame("Missing Particle");
		frameX.setSize(1000, 1000);
		EmbeddedCanvas canX = new EmbeddedCanvas();
		frameX.add(canX);
		frameX.setLocationRelativeTo(null);
		frameX.setVisible(true);
		canX.divide(2, 2);
		canX.cd(0);
		canX.setFont("Arial");
		htheta_cone_gamma_rb.setTitle("#Delta#theta_cone(#gamma)");
		htheta_cone_gamma_rb.setTitleX("#Delta#theta_cone(#gamma) [deg]");
		htheta_cone_gamma_rb.setTitleY("Counts");
	//	htheta_cone_gamma_rb.setLineWidth(2);
		htheta_cone_gamma_rb.setOptStat(10);
		canX.getPad(0).setTitleFontSize(32);
		canX.getPad(0).setAxisTitleFontSize(32);
		canX.getPad(0).setAxisLabelFontSize(24);
		canX.getPad(0).setStatBoxFontSize(18);
		canX.draw(htheta_cone_gamma_rb, "same");
		canX.cd(1);
		canX.setFont("Arial");
		hX_proton_m_rb.setTitle("M_X_(ep->e'#gamma)");
		hX_proton_m_rb.setTitleX("M_X [GeV]");
		hX_proton_m_rb.setTitleY("Counts");
	//	hX_proton_m_rb.setLineWidth(2);
		hX_proton_m_rb.setOptStat(10);
		canX.getPad(1).setTitleFontSize(32);
		canX.getPad(1).setAxisTitleFontSize(32);
		canX.getPad(1).setAxisLabelFontSize(24);
		canX.getPad(1).setStatBoxFontSize(18);
		canX.draw(hX_proton_m_rb, "same");
		canX.cd(2);
		canX.setFont("Arial");
		hX_epg_pt_rb.setTitle("p_X_(ep->e'p'#gamma)");
		hX_epg_pt_rb.setTitleX("p_X [GeV]");
		hX_epg_pt_rb.setTitleY("Counts");
	//	hX_epg_pt_rb.setLineWidth(2);
		hX_epg_pt_rb.setOptStat(10);
		canX.getPad(2).setTitleFontSize(32);
		canX.getPad(2).setAxisTitleFontSize(32);
		canX.getPad(2).setAxisLabelFontSize(24);
		canX.getPad(2).setStatBoxFontSize(18);
		canX.draw(hX_epg_pt_rb, "same");
		canX.cd(3);
		canX.setFont("Arial");
		hX_epg_E_rb.setTitle("E_X_(ep->e'p'#gamma)");
		hX_epg_E_rb.setTitleX("E_X [GeV]");
		hX_epg_E_rb.setTitleY("Counts");
	//	hX_epg_E_rb.setLineWidth(2);
		hX_epg_E_rb.setOptStat(10);
		canX.getPad(3).setTitleFontSize(32);
		canX.getPad(3).setAxisTitleFontSize(32);
		canX.getPad(3).setAxisLabelFontSize(24);
		canX.getPad(3).setStatBoxFontSize(18);
		canX.draw(hX_epg_E_rb, "same");
		
		JFrame framekin = new JFrame("Kinematic Variables");
		framekin.setSize(1000, 500);
		EmbeddedCanvas cankin = new EmbeddedCanvas();
		framekin.add(cankin);
		framekin.setLocationRelativeTo(null);
		framekin.setVisible(true);
		cankin.divide(2, 1);
		cankin.cd(0);
		cankin.setFont("Arial");
		hQ2_vs_x_B_rb.setTitle("Q^2 vs. x_B");
		hQ2_vs_x_B_rb.setTitleX("x_B");
		hQ2_vs_x_B_rb.setTitleY("Q^2 [GeV^2]");
		cankin.getPad(0).setTitleFontSize(32);
		cankin.getPad(0).setAxisTitleFontSize(32);
		cankin.getPad(0).setAxisLabelFontSize(24);
		cankin.getPad(0).setStatBoxFontSize(18);
		cankin.draw(hQ2_vs_x_B_rb, "same");
		cankin.cd(1);
		cankin.setFont("Arial");
		hnegt_vs_phi_trento.setTitle("-t vs. #phi");
		hnegt_vs_phi_trento.setTitleX("#phi [deg]");
		hnegt_vs_phi_trento.setTitleY("-t [GeV^2]");
		cankin.getPad(1).setTitleFontSize(32);
		cankin.getPad(1).setAxisTitleFontSize(32);
		cankin.getPad(1).setAxisLabelFontSize(24);
		cankin.getPad(1).setStatBoxFontSize(18);
		cankin.draw(hnegt_vs_phi_trento, "same");
		
		JFrame framebsa = new JFrame("BSA");
		framebsa.setSize(1500, 500);
		EmbeddedCanvas canbsa = new EmbeddedCanvas();
		framebsa.add(canbsa);
		framebsa.setLocationRelativeTo(null);
		framebsa.setVisible(true);
		canbsa.divide(3, 1);
		canbsa.cd(0);
		canbsa.setFont("Arial");
		hphi_CM_poshel.setTitle("(+) Helicity #phi");
		hphi_CM_poshel.setTitleX("#phi [deg]");
		hphi_CM_poshel.setTitleY("Counts");
	//	hphi_CM_poshel.setLineWidth(2);
		hphi_CM_poshel.setOptStat(10);
		canbsa.getPad(0).setTitleFontSize(32);
		canbsa.getPad(0).setAxisTitleFontSize(32);
		canbsa.getPad(0).setAxisLabelFontSize(24);
		canbsa.getPad(0).setStatBoxFontSize(18);
		canbsa.draw(hphi_CM_poshel, "same");
		canbsa.cd(1);
		canbsa.setFont("Arial");
		hphi_CM_neghel.setTitle("(-) Helicity #phi");
		hphi_CM_neghel.setTitleX("#phi [deg]");
		hphi_CM_neghel.setTitleY("Counts");
	//	hphi_CM_neghel.setLineWidth(2);
		hphi_CM_neghel.setOptStat(10);
		canbsa.getPad(1).setTitleFontSize(32);
		canbsa.getPad(1).setAxisTitleFontSize(32);
		canbsa.getPad(1).setAxisLabelFontSize(24);
		canbsa.getPad(1).setStatBoxFontSize(18);
		canbsa.draw(hphi_CM_neghel, "same");
		canbsa.cd(2);
		canbsa.setFont("Arial");
		hBSAy.setTitle("Raw Beam Spin Asymmetry");
		hBSAy.setTitleX("#phi [deg]");
		hBSAy.setTitleY("Raw BSA");
		hBSAy.setMarkerSize(5);
		hBSAy.setMarkerStyle(1);
		hBSAy.setLineThickness(1);
		canbsa.getPad(2).setTitleFontSize(32);
		canbsa.getPad(2).setAxisTitleFontSize(32);
		canbsa.getPad(2).setAxisLabelFontSize(24);
		canbsa.getPad(2).setAxisRange(0, 360, -0.3, 0.3);
		canbsa.getPad(2).setStatBoxFontSize(18);
		canbsa.draw(hBSAy, "same");
		F1D fBSA= new F1D("fBSA", "[A]*sin(x/57.3)/(1+[B]*cos(x/57.3))", 0, 360);
		fBSA.setParameter(0, 0.1);
		fBSA.setParameter(1, -0.5);
		DataFitter.fit(fBSA, hBSAy, "Q");
		fBSA.setLineColor(2);
		fBSA.setLineWidth(3);
		fBSA.setOptStat(1110);
		canbsa.draw(fBSA, "same");
		System.out.println("==========================================================================================");
		System.out.println("BSA:	x	#deltax	y	#delta	y");
		for(int p = 0; p < 30; p++)		
		{
			System.out.println("	" + phibc[p] + "	" + 0 + "	" + bsa.get(p) + "	" + bsaerr.get(p));
		}
		System.out.println("A: " + fBSA.parameter(0).value() + " +/- " + fBSA.parameter(0).error());
		System.out.println("B: " + fBSA.parameter(1).value() + " +/- " + fBSA.parameter(1).error());
		
		JFrame frame_phi_poshel_bin = new JFrame("Binned (+) Helicity #phi Distribution");
		frame_phi_poshel_bin.setSize(1200, 600);
		EmbeddedCanvas can_phi_poshel_bin = new EmbeddedCanvas();
		frame_phi_poshel_bin.add(can_phi_poshel_bin);
		frame_phi_poshel_bin.setLocationRelativeTo(null);
		frame_phi_poshel_bin.setVisible(true);
		can_phi_poshel_bin.divide(4, 2);
		for(int ih = 0; ih  < 8; ih++)
		{
			can_phi_poshel_bin.cd(ih);
			can_phi_poshel_bin.setFont("Arial");
			histGroups_phi_CM_poshel_bin.getItem(ih).setTitle("Bin " + (ih+1));
			histGroups_phi_CM_poshel_bin.getItem(ih).setTitleX("#phi [#degree]");
			histGroups_phi_CM_poshel_bin.getItem(ih).setOptStat(10);
			can_phi_poshel_bin.getPad(ih).setAxisLabelFontSize(10);
			can_phi_poshel_bin.getPad(ih).setAxisTitleFontSize(10);
			can_phi_poshel_bin.getPad(ih).setTitleFontSize(10);
			can_phi_poshel_bin.draw(histGroups_phi_CM_poshel_bin.getItem(ih), "same");
		}
		
		JFrame frame_phi_neghel_bin = new JFrame("Binned (-) Helicity #phi Distribution");
		frame_phi_neghel_bin.setSize(1200, 600);
		EmbeddedCanvas can_phi_neghel_bin = new EmbeddedCanvas();
		frame_phi_neghel_bin.add(can_phi_neghel_bin);
		frame_phi_neghel_bin.setLocationRelativeTo(null);
		frame_phi_neghel_bin.setVisible(true);
		can_phi_neghel_bin.divide(4, 2);
		for(int ih = 0; ih  < 8; ih++)
		{
			can_phi_neghel_bin.cd(ih);
			can_phi_neghel_bin.setFont("Arial");
			histGroups_phi_CM_neghel_bin.getItem(ih).setTitle("Bin " + (ih+1));
			histGroups_phi_CM_neghel_bin.getItem(ih).setTitleX("#phi [#degree]");
			histGroups_phi_CM_neghel_bin.getItem(ih).setOptStat(10);
			can_phi_neghel_bin.getPad(ih).setAxisLabelFontSize(10);
			can_phi_neghel_bin.getPad(ih).setAxisTitleFontSize(10);
			can_phi_neghel_bin.getPad(ih).setTitleFontSize(10);
			can_phi_neghel_bin.draw(histGroups_phi_CM_neghel_bin.getItem(ih), "same");
		}
		
		JFrame framebsa_bin = new JFrame("Binned BSA");
		framebsa_bin.setSize(1200, 600);
		EmbeddedCanvas canbsa_bin = new EmbeddedCanvas();
		framebsa_bin.add(canbsa_bin);
		framebsa_bin.setLocationRelativeTo(null);
		framebsa_bin.setVisible(true);
		canbsa_bin.divide(4, 2);
		for(int ih = 0; ih  < 8; ih++)
		{
			canbsa_bin.cd(ih);
			canbsa_bin.setFont("Arial");
			hBSAy_bin.getItem(ih).setMarkerSize(3);
			hBSAy_bin.getItem(ih).setMarkerStyle(1);
			hBSAy_bin.getItem(ih).setLineThickness(1);
			hBSAy_bin.getItem(ih).setTitle("Bin " + (ih+1));
			hBSAy_bin.getItem(ih).setTitleX("#phi [#degree]");
			hBSAy_bin.getItem(ih).setTitleY("Raw BSA");
			canbsa_bin.getPad(ih).setAxisLabelFontSize(10);
			canbsa_bin.getPad(ih).setAxisTitleFontSize(10);
			canbsa_bin.getPad(ih).setTitleFontSize(10);
			canbsa_bin.getPad(ih).setAxisRange(0, 360, -0.3, 0.3);
			canbsa_bin.draw(hBSAy_bin.getItem(ih), "same");
			F1D fBSA_bin= new F1D("fBSA_bin", "[A]*sin(x/57.3)/(1+[B]*cos(x/57.3))", 0, 360);
			fBSA_bin.setParameter(0, 0.15);
			fBSA_bin.setParameter(1, -0.5);
			DataFitter.fit(fBSA_bin, hBSAy_bin.getItem(ih), "Q");
			fBSA_bin.setLineColor(2);
			fBSA_bin.setLineWidth(2);
			fBSA_bin.setOptStat(1110);
			canbsa_bin.draw(fBSA_bin, "same");
			System.out.println("==========================================================================================");
			System.out.println("BSA_bin:	x	#deltax	y	#delta	y");
			System.out.println("Bin No.: " + ih);
			for(int p = 0; p < 30; p++)		
			{
				System.out.println("	" + phibc[p] + "	" + 0 + "	" + bsa_bin.getItem(ih, p) + "	" + bsaerr_bin.getItem(ih, p));
			}
			System.out.println("A: " + fBSA_bin.parameter(0).value() + " +/- " + fBSA_bin.parameter(0).error());
			System.out.println("B: " + fBSA_bin.parameter(1).value() + " +/- " + fBSA_bin.parameter(1).error());
	//		System.out.println("A: " + fBSA_bin.parameter(0).value() + " +/- " + fBSA_bin.parameter(0).error());
	//		System.out.println("B: " + fBSA_bin.parameter(1).value() + " +/- " + fBSA_bin.parameter(1).error());
		}
		
		JFrame frame_alphat = new JFrame("#alpha vs. -t");
		frame_alphat.setSize(1000, 1000);
		EmbeddedCanvas can_alphat = new EmbeddedCanvas();
		frame_alphat.add(can_alphat);
		frame_alphat.setLocationRelativeTo(null);
		frame_alphat.setVisible(true);
		can_alphat.divide(2, 2);
		can_alphat.cd(0);
		can_alphat.setFont("Arial");
		halpha_vs_t.getItem(0).setTitle("Raw BSA #alpha vs. -t");
		halpha_vs_t.getItem(0).setTitleX("-t [GeV^2]");
		halpha_vs_t.getItem(0).setTitleY("Pre-Correction #alpha");
		halpha_vs_t.getItem(0).setMarkerSize(5);
		halpha_vs_t.getItem(0).setMarkerStyle(1);
		halpha_vs_t.getItem(0).setLineThickness(1);
		can_alphat.getPad(0).setTitleFontSize(32);
		can_alphat.getPad(0).setAxisTitleFontSize(32);
		can_alphat.getPad(0).setAxisLabelFontSize(24);
		can_alphat.getPad(0).setAxisRange(0, 2, -0.1, 0.45);
		can_alphat.getPad(0).setStatBoxFontSize(18);
		can_alphat.draw(halpha_vs_t.getItem(0), "same");
		can_alphat.cd(1);
		can_alphat.setFont("Arial");
		halpha_vs_t.getItem(1).setTitle("Raw BSA #alpha vs. -t");
		halpha_vs_t.getItem(1).setTitleX("-t [GeV^2]");
		halpha_vs_t.getItem(1).setTitleY("Pre-Correction #alpha");
		halpha_vs_t.getItem(1).setMarkerSize(5);
		halpha_vs_t.getItem(1).setMarkerStyle(1);
		halpha_vs_t.getItem(1).setLineThickness(1);
		can_alphat.getPad(1).setTitleFontSize(32);
		can_alphat.getPad(1).setAxisTitleFontSize(32);
		can_alphat.getPad(1).setAxisLabelFontSize(24);
		can_alphat.getPad(1).setAxisRange(0, 2, -0.1, 0.45);
		can_alphat.getPad(1).setStatBoxFontSize(18);
		can_alphat.draw(halpha_vs_t.getItem(1), "same");
		halpha_vs_t.getItem(2).setMarkerSize(5);
		halpha_vs_t.getItem(2).setLineThickness(1);
		halpha_vs_t.getItem(2).setMarkerStyle(0);
		halpha_vs_t.getItem(2).setMarkerColor(2);
		halpha_vs_t.getItem(2).setLineColor(2);
		can_alphat.draw(halpha_vs_t.getItem(2), "same");
		can_alphat.cd(2);
		can_alphat.setFont("Arial");
		halpha_vs_t.getItem(3).setTitle("Raw BSA #alpha vs. -t");
		halpha_vs_t.getItem(3).setTitleX("-t [GeV^2]");
		halpha_vs_t.getItem(3).setTitleY("Pre-Correction #alpha");
		halpha_vs_t.getItem(3).setMarkerSize(5);
		halpha_vs_t.getItem(3).setMarkerStyle(1);
		halpha_vs_t.getItem(3).setLineThickness(1);
		can_alphat.getPad(2).setTitleFontSize(32);
		can_alphat.getPad(2).setAxisTitleFontSize(32);
		can_alphat.getPad(2).setAxisLabelFontSize(24);
		can_alphat.getPad(2).setAxisRange(0, 2, -0.1, 0.45);
		can_alphat.getPad(2).setStatBoxFontSize(18);
		can_alphat.draw(halpha_vs_t.getItem(3), "same");
		halpha_vs_t.getItem(4).setMarkerSize(5);
		halpha_vs_t.getItem(4).setLineThickness(1);
		halpha_vs_t.getItem(4).setMarkerStyle(0);
		halpha_vs_t.getItem(4).setMarkerColor(2);
		halpha_vs_t.getItem(4).setLineColor(2);
		can_alphat.draw(halpha_vs_t.getItem(4), "same");
		halpha_vs_t.getItem(5).setMarkerSize(5);
		halpha_vs_t.getItem(5).setLineThickness(1);
		halpha_vs_t.getItem(5).setMarkerStyle(2);
		halpha_vs_t.getItem(5).setMarkerColor(4);
		halpha_vs_t.getItem(5).setLineColor(4);
		can_alphat.draw(halpha_vs_t.getItem(5), "same");
		can_alphat.cd(3);
		can_alphat.setFont("Arial");
		halpha_vs_t.getItem(6).setTitle("Raw BSA #alpha vs. -t'");
		halpha_vs_t.getItem(6).setTitleX("-t [GeV^2]");
		halpha_vs_t.getItem(6).setTitleY("Pre-Correction #alpha");
		halpha_vs_t.getItem(6).setMarkerSize(5);
		halpha_vs_t.getItem(6).setMarkerStyle(1);
		halpha_vs_t.getItem(6).setLineThickness(1);
		can_alphat.getPad(3).setTitleFontSize(32);
		can_alphat.getPad(3).setAxisTitleFontSize(32);
		can_alphat.getPad(3).setAxisLabelFontSize(24);
		can_alphat.getPad(3).setAxisRange(0, 2, -0.1, 0.45);
		can_alphat.getPad(3).setStatBoxFontSize(18);
		can_alphat.draw(halpha_vs_t.getItem(6), "same");
		halpha_vs_t.getItem(7).setMarkerSize(5);
		halpha_vs_t.getItem(7).setLineThickness(1);
		halpha_vs_t.getItem(7).setMarkerStyle(0);
		halpha_vs_t.getItem(7).setMarkerColor(2);
		halpha_vs_t.getItem(7).setLineColor(2);
		can_alphat.draw(halpha_vs_t.getItem(7), "same");
		
		JFrame framepi0X = new JFrame("Missing #pi^0");
		framepi0X.setSize(1000, 1000);
		EmbeddedCanvas canpi0X = new EmbeddedCanvas();
		framepi0X.add(canpi0X);
		framepi0X.setLocationRelativeTo(null);
		framepi0X.setVisible(true);
		canpi0X.divide(2, 2);
		canpi0X.cd(0);
		canpi0X.setFont("Arial");
		htheta_cone_pi0_rb.setTitle("#Delta#theta_cone(#pi^0)");
		htheta_cone_pi0_rb.setTitleX("#Delta#theta_cone(#pi^0) [deg]");
		htheta_cone_pi0_rb.setTitleY("Counts");
	//	htheta_cone_pi0_rb.setLineWidth(2);
		htheta_cone_pi0_rb.setOptStat(10);
		canpi0X.getPad(0).setTitleFontSize(32);
		canpi0X.getPad(0).setAxisTitleFontSize(32);
		canpi0X.getPad(0).setAxisLabelFontSize(24);
		canpi0X.getPad(0).setStatBoxFontSize(18);
		canpi0X.draw(htheta_cone_pi0_rb, "same");
		canpi0X.cd(1);
		canpi0X.setFont("Arial");
		hX_proton_m_pi0_rb.setTitle("M_X_(ep->e'#gamma#gamma)");
		hX_proton_m_pi0_rb.setTitleX("M_X [GeV]");
		hX_proton_m_pi0_rb.setTitleY("Counts");
	//	hX_proton_m_pi0_rb.setLineWidth(2);
		hX_proton_m_pi0_rb.setOptStat(10);
		canpi0X.getPad(1).setTitleFontSize(32);
		canpi0X.getPad(1).setAxisTitleFontSize(32);
		canpi0X.getPad(1).setAxisLabelFontSize(24);
		canpi0X.getPad(1).setStatBoxFontSize(18);
		canpi0X.draw(hX_proton_m_pi0_rb, "same");
		canpi0X.cd(2);
		canpi0X.setFont("Arial");
		hX_epgg_pt_rb.setTitle("p_X_(ep->e'p'#gamma#gamma)");
		hX_epgg_pt_rb.setTitleX("p_X [GeV]");
		hX_epgg_pt_rb.setTitleY("Counts");
	//	hX_epgg_pt_rb.setLineWidth(2);
		hX_epgg_pt_rb.setOptStat(10);
		canpi0X.getPad(2).setTitleFontSize(32);
		canpi0X.getPad(2).setAxisTitleFontSize(32);
		canpi0X.getPad(2).setAxisLabelFontSize(24);
		canpi0X.getPad(2).setStatBoxFontSize(18);
		canpi0X.draw(hX_epgg_pt_rb, "same");
		canpi0X.cd(3);
		canpi0X.setFont("Arial");
		hX_epgg_E_rb.setTitle("E_X_(ep->e'p'#gamma#gamma)");
		hX_epgg_E_rb.setTitleX("E_X [GeV]");
		hX_epgg_E_rb.setTitleY("Counts");
	//	hX_epgg_E_rb.setLineWidth(2);
		hX_epgg_E_rb.setOptStat(10);
		canpi0X.getPad(3).setTitleFontSize(32);
		canpi0X.getPad(3).setAxisTitleFontSize(32);
		canpi0X.getPad(3).setAxisLabelFontSize(24);
		canpi0X.getPad(3).setStatBoxFontSize(18);
		canpi0X.draw(hX_epgg_E_rb, "same");
		
		JFrame framebsa_pi0 = new JFrame("#pi^0 BSA");
		framebsa_pi0.setSize(1500, 500);
		EmbeddedCanvas canbsa_pi0 = new EmbeddedCanvas();
		framebsa_pi0.add(canbsa_pi0);
		framebsa_pi0.setLocationRelativeTo(null);
		framebsa_pi0.setVisible(true);
		canbsa_pi0.divide(3, 1);
		canbsa_pi0.cd(0);
		canbsa_pi0.setFont("Arial");
		hphi_CM_poshel_pi0.setTitle("(+) Helicity #phi");
		hphi_CM_poshel_pi0.setTitleX("#phi [deg]");
		hphi_CM_poshel_pi0.setTitleY("Counts");
	//	hphi_CM_poshel_pi0.setLineWidth(2);
		hphi_CM_poshel_pi0.setOptStat(10);
		canbsa_pi0.getPad(0).setTitleFontSize(32);
		canbsa_pi0.getPad(0).setAxisTitleFontSize(32);
		canbsa_pi0.getPad(0).setAxisLabelFontSize(24);
		canbsa_pi0.getPad(0).setStatBoxFontSize(18);
		canbsa_pi0.draw(hphi_CM_poshel_pi0, "same");
		canbsa_pi0.cd(1);
		canbsa_pi0.setFont("Arial");
		hphi_CM_neghel_pi0.setTitle("(-) Helicity #phi");
		hphi_CM_neghel_pi0.setTitleX("#phi [deg]");
		hphi_CM_neghel_pi0.setTitleY("Counts");
	//	hphi_CM_neghel_pi0.setLineWidth(2);
		hphi_CM_neghel_pi0.setOptStat(10);
		canbsa_pi0.getPad(1).setTitleFontSize(32);
		canbsa_pi0.getPad(1).setAxisTitleFontSize(32);
		canbsa_pi0.getPad(1).setAxisLabelFontSize(24);
		canbsa_pi0.getPad(1).setStatBoxFontSize(18);
		canbsa_pi0.draw(hphi_CM_neghel_pi0, "same");
		canbsa_pi0.cd(2);
		canbsa_pi0.setFont("Arial");
		pi0_hBSAy.setTitle("#pi^0 Beam Spin Asymmetry");
		pi0_hBSAy.setTitleX("#phi [deg]");
		pi0_hBSAy.setTitleY("#pi^0 BSA");
		pi0_hBSAy.setMarkerSize(5);
		pi0_hBSAy.setMarkerStyle(1);
		pi0_hBSAy.setLineThickness(1);
		canbsa_pi0.getPad(2).setTitleFontSize(32);
		canbsa_pi0.getPad(2).setAxisTitleFontSize(32);
		canbsa_pi0.getPad(2).setAxisLabelFontSize(24);
		canbsa_pi0.getPad(2).setAxisRange(0, 360, -0.3, 0.3);
		canbsa_pi0.getPad(2).setStatBoxFontSize(18);
		canbsa_pi0.draw(pi0_hBSAy, "same");
		F1D fBSA_pi0= new F1D("fBSA_pi0", "[A]*sin(x/57.3)", 0, 360);
		fBSA_pi0.setParameter(0, 0.05);
		DataFitter.fit(fBSA_pi0, pi0_hBSAy, "Q");
		fBSA_pi0.setLineColor(2);
		fBSA_pi0.setLineWidth(3);
		fBSA_pi0.setOptStat(110);
		canbsa_pi0.draw(fBSA_pi0, "same");
		System.out.println("==========================================================================================");
		System.out.println("BSA_pi0:	x	#deltax	y	#delta	y");
		for(int p = 0; p < 30; p++)		
		{
			System.out.println("	" + phibc[p] + "	" + 0 + "	" + pi0_bsa.get(p) + "	" + pi0_bsaerr.get(p));
		}
		System.out.println("A: " + fBSA_pi0.parameter(0).value() + " +/- " + fBSA_pi0.parameter(0).error());
		
		JFrame frame_phi_poshel_pi0_bin = new JFrame("Binned (+) Helicity #pi^0 #phi Distribution");
		frame_phi_poshel_pi0_bin.setSize(1200, 600);
		EmbeddedCanvas can_phi_poshel_pi0_bin = new EmbeddedCanvas();
		frame_phi_poshel_pi0_bin.add(can_phi_poshel_pi0_bin);
		frame_phi_poshel_pi0_bin.setLocationRelativeTo(null);
		frame_phi_poshel_pi0_bin.setVisible(true);
		can_phi_poshel_pi0_bin.divide(4, 2);
		for(int ih = 0; ih  < 8; ih++)
		{
			can_phi_poshel_pi0_bin.cd(ih);
			can_phi_poshel_pi0_bin.setFont("Arial");
			histGroups_phi_CM_poshel_pi0_bin.getItem(ih).setTitle("Bin " + (ih+1));
			histGroups_phi_CM_poshel_pi0_bin.getItem(ih).setTitleX("#phi [#degree]");
			histGroups_phi_CM_poshel_pi0_bin.getItem(ih).setOptStat(10);
			can_phi_poshel_pi0_bin.getPad(ih).setAxisLabelFontSize(10);
			can_phi_poshel_pi0_bin.getPad(ih).setAxisTitleFontSize(10);
			can_phi_poshel_pi0_bin.getPad(ih).setTitleFontSize(10);
			can_phi_poshel_pi0_bin.draw(histGroups_phi_CM_poshel_pi0_bin.getItem(ih), "same");
		}
		
		JFrame frame_phi_neghel_pi0_bin = new JFrame("Binned (-) Helicity #pi^0 #phi Distribution");
		frame_phi_neghel_pi0_bin.setSize(1200, 600);
		EmbeddedCanvas can_phi_neghel_pi0_bin = new EmbeddedCanvas();
		frame_phi_neghel_pi0_bin.add(can_phi_neghel_pi0_bin);
		frame_phi_neghel_pi0_bin.setLocationRelativeTo(null);
		frame_phi_neghel_pi0_bin.setVisible(true);
		can_phi_neghel_pi0_bin.divide(4, 2);
		for(int ih = 0; ih  < 8; ih++)
		{
			can_phi_neghel_pi0_bin.cd(ih);
			can_phi_neghel_pi0_bin.setFont("Arial");
			histGroups_phi_CM_neghel_pi0_bin.getItem(ih).setTitle("Bin " + (ih+1));
			histGroups_phi_CM_neghel_pi0_bin.getItem(ih).setTitleX("#phi [#degree]");
			histGroups_phi_CM_neghel_pi0_bin.getItem(ih).setOptStat(10);
			can_phi_neghel_pi0_bin.getPad(ih).setAxisLabelFontSize(10);
			can_phi_neghel_pi0_bin.getPad(ih).setAxisTitleFontSize(10);
			can_phi_neghel_pi0_bin.getPad(ih).setTitleFontSize(10);
			can_phi_neghel_pi0_bin.draw(histGroups_phi_CM_neghel_pi0_bin.getItem(ih), "same");
		}
		
		JFrame framebsa_pi0_bin = new JFrame("Binned #pi^0 BSA");
		framebsa_pi0_bin.setSize(1200, 600);
		EmbeddedCanvas canbsa_pi0_bin = new EmbeddedCanvas();
		framebsa_pi0_bin.add(canbsa_pi0_bin);
		framebsa_pi0_bin.setLocationRelativeTo(null);
		framebsa_pi0_bin.setVisible(true);
		canbsa_pi0_bin.divide(4, 2);
		for(int ih = 0; ih  < 8; ih++)
		{
			canbsa_pi0_bin.cd(ih);
			canbsa_pi0_bin.setFont("Arial");
			pi0_hBSAy_bin.getItem(ih).setMarkerSize(3);
			pi0_hBSAy_bin.getItem(ih).setMarkerStyle(1);
			pi0_hBSAy_bin.getItem(ih).setLineThickness(1);
			pi0_hBSAy_bin.getItem(ih).setTitle("Bin " + (ih+1));
			pi0_hBSAy_bin.getItem(ih).setTitleX("#phi [#degree]");
			pi0_hBSAy_bin.getItem(ih).setTitleY("#pi^0 BSA");
			canbsa_pi0_bin.getPad(ih).setAxisLabelFontSize(10);
			canbsa_pi0_bin.getPad(ih).setAxisTitleFontSize(10);
			canbsa_pi0_bin.getPad(ih).setTitleFontSize(10);
			canbsa_pi0_bin.getPad(ih).setAxisRange(0, 360, -0.3, 0.3);
			canbsa_pi0_bin.draw(pi0_hBSAy_bin.getItem(ih), "same");
			F1D fBSA_bin= new F1D("fBSA_bin", "[A]*sin(x/57.3)", 0, 360);
			fBSA_bin.setParameter(0, 0.05);
			DataFitter.fit(fBSA_bin, pi0_hBSAy_bin.getItem(ih), "Q");
			fBSA_bin.setLineColor(2);
			fBSA_bin.setOptStat(110);
			canbsa_pi0_bin.draw(fBSA_bin, "same");
	//		System.out.println("Bin No.: " + ih);
	//		System.out.println("A: " + fBSA_bin.parameter(0).value() + " +/- " + fBSA_bin.parameter(0).error());
		}
		
		JFrame frame_alphat_pi0 = new JFrame("#pi^0 #alpha vs. -t");
		frame_alphat_pi0.setSize(1000, 1000);
		EmbeddedCanvas can_alphat_pi0 = new EmbeddedCanvas();
		frame_alphat_pi0.add(can_alphat_pi0);
		frame_alphat_pi0.setLocationRelativeTo(null);
		frame_alphat_pi0.setVisible(true);
		can_alphat_pi0.divide(2, 2);
		can_alphat_pi0.cd(0);
		can_alphat_pi0.setFont("Arial");
		pi0_halpha_vs_t.getItem(0).setTitle("#pi^0 BSA #alpha vs. -t");
		pi0_halpha_vs_t.getItem(0).setTitleX("-t [GeV^2]");
		pi0_halpha_vs_t.getItem(0).setTitleY("#pi^0 #alpha");
		pi0_halpha_vs_t.getItem(0).setMarkerSize(5);
		pi0_halpha_vs_t.getItem(0).setMarkerStyle(1);
		pi0_halpha_vs_t.getItem(0).setLineThickness(1);
		can_alphat_pi0.getPad(0).setTitleFontSize(32);
		can_alphat_pi0.getPad(0).setAxisTitleFontSize(32);
		can_alphat_pi0.getPad(0).setAxisLabelFontSize(24);
		can_alphat_pi0.getPad(0).setAxisRange(0, 2, -0.1, 0.45);
		can_alphat_pi0.getPad(0).setStatBoxFontSize(18);
		can_alphat_pi0.draw(pi0_halpha_vs_t.getItem(0), "same");
		can_alphat_pi0.cd(1);
		can_alphat_pi0.setFont("Arial");
		pi0_halpha_vs_t.getItem(1).setTitle("#pi^0 BSA #alpha vs. -t");
		pi0_halpha_vs_t.getItem(1).setTitleX("-t [GeV^2]");
		pi0_halpha_vs_t.getItem(1).setTitleY("#pi^0 #alpha");
		pi0_halpha_vs_t.getItem(1).setMarkerSize(5);
		pi0_halpha_vs_t.getItem(1).setMarkerStyle(1);
		pi0_halpha_vs_t.getItem(1).setLineThickness(1);
		can_alphat_pi0.getPad(1).setTitleFontSize(32);
		can_alphat_pi0.getPad(1).setAxisTitleFontSize(32);
		can_alphat_pi0.getPad(1).setAxisLabelFontSize(24);
		can_alphat_pi0.getPad(1).setAxisRange(0, 2, -0.1, 0.45);
		can_alphat_pi0.getPad(1).setStatBoxFontSize(18);
		can_alphat_pi0.draw(pi0_halpha_vs_t.getItem(1), "same");
		pi0_halpha_vs_t.getItem(2).setMarkerSize(5);
		pi0_halpha_vs_t.getItem(2).setLineThickness(1);
		pi0_halpha_vs_t.getItem(2).setMarkerStyle(0);
		pi0_halpha_vs_t.getItem(2).setMarkerColor(2);
		pi0_halpha_vs_t.getItem(2).setLineColor(2);
		can_alphat_pi0.draw(pi0_halpha_vs_t.getItem(2), "same");
		can_alphat_pi0.cd(2);
		can_alphat_pi0.setFont("Arial");
		pi0_halpha_vs_t.getItem(3).setTitle("#pi^0 BSA #alpha vs. -t");
		pi0_halpha_vs_t.getItem(3).setTitleX("-t [GeV^2]");
		pi0_halpha_vs_t.getItem(3).setTitleY("#pi^0 #alpha");
		pi0_halpha_vs_t.getItem(3).setMarkerSize(5);
		pi0_halpha_vs_t.getItem(3).setMarkerStyle(1);
		pi0_halpha_vs_t.getItem(3).setLineThickness(1);
		can_alphat_pi0.getPad(2).setTitleFontSize(32);
		can_alphat_pi0.getPad(2).setAxisTitleFontSize(32);
		can_alphat_pi0.getPad(2).setAxisLabelFontSize(24);
		can_alphat_pi0.getPad(2).setAxisRange(0, 2, -0.1, 0.45);
		can_alphat_pi0.getPad(2).setStatBoxFontSize(18);
		can_alphat_pi0.draw(pi0_halpha_vs_t.getItem(3), "same");
		pi0_halpha_vs_t.getItem(4).setMarkerSize(5);
		pi0_halpha_vs_t.getItem(4).setLineThickness(1);
		pi0_halpha_vs_t.getItem(4).setMarkerStyle(0);
		pi0_halpha_vs_t.getItem(4).setMarkerColor(2);
		pi0_halpha_vs_t.getItem(4).setLineColor(2);
		can_alphat_pi0.draw(pi0_halpha_vs_t.getItem(4), "same");
		pi0_halpha_vs_t.getItem(5).setMarkerSize(5);
		pi0_halpha_vs_t.getItem(5).setLineThickness(1);
		pi0_halpha_vs_t.getItem(5).setMarkerStyle(2);
		pi0_halpha_vs_t.getItem(5).setMarkerColor(4);
		pi0_halpha_vs_t.getItem(5).setLineColor(4);
		can_alphat_pi0.draw(pi0_halpha_vs_t.getItem(5), "same");
		can_alphat_pi0.cd(3);
		can_alphat_pi0.setFont("Arial");
		pi0_halpha_vs_t.getItem(6).setTitle("#pi^0 BSA #alpha vs. -t");
		pi0_halpha_vs_t.getItem(6).setTitleX("-t [GeV^2]");
		pi0_halpha_vs_t.getItem(6).setTitleY("#pi^0 #alpha");
		pi0_halpha_vs_t.getItem(6).setMarkerSize(5);
		pi0_halpha_vs_t.getItem(6).setMarkerStyle(1);
		pi0_halpha_vs_t.getItem(6).setLineThickness(1);
		can_alphat_pi0.getPad(3).setTitleFontSize(32);
		can_alphat_pi0.getPad(3).setAxisTitleFontSize(32);
		can_alphat_pi0.getPad(3).setAxisLabelFontSize(24);
		can_alphat_pi0.getPad(3).setAxisRange(0, 2, -0.1, 0.45);
		can_alphat_pi0.getPad(3).setStatBoxFontSize(18);
		can_alphat_pi0.draw(pi0_halpha_vs_t.getItem(6), "same");
		pi0_halpha_vs_t.getItem(7).setMarkerSize(5);
		pi0_halpha_vs_t.getItem(7).setLineThickness(1);
		pi0_halpha_vs_t.getItem(7).setMarkerStyle(0);
		pi0_halpha_vs_t.getItem(7).setMarkerColor(2);
		pi0_halpha_vs_t.getItem(7).setLineColor(2);
		can_alphat_pi0.draw(pi0_halpha_vs_t.getItem(7), "same");
		
		
		JFrame framer_pi0_sim = new JFrame("#pi^0_sim r");
		framer_pi0_sim.setSize(500, 500);
		EmbeddedCanvas canr_pi0_sim = new EmbeddedCanvas();
		framer_pi0_sim.add(canr_pi0_sim);
		framer_pi0_sim.setLocationRelativeTo(null);
		framer_pi0_sim.setVisible(true);
		canr_pi0_sim.divide(1, 1);
		canr_pi0_sim.cd(0);
		canr_pi0_sim.setFont("Arial");
		pi0_sim_hr.setTitle("r vs. #phi");
		pi0_sim_hr.setTitleX("#phi [deg]");
		pi0_sim_hr.setTitleY("r");
	//	pi0_sim_hr.setLineWidth(2);
		pi0_sim_hr.setMarkerSize(5);
		pi0_sim_hr.setMarkerStyle(1);
		canr_pi0_sim.getPad(0).setTitleFontSize(32);
		canr_pi0_sim.getPad(0).setAxisTitleFontSize(32);
		canr_pi0_sim.getPad(0).setAxisLabelFontSize(24);
		canr_pi0_sim.getPad(0).setAxisRange(0, 360, 0, 2);
		canr_pi0_sim.draw(pi0_sim_hr, "same");
		
		JFrame framef_pi0_sim = new JFrame("#pi^0 f");
		framef_pi0_sim.setSize(500, 500);
		EmbeddedCanvas canf_pi0_sim = new EmbeddedCanvas();
		framef_pi0_sim.add(canf_pi0_sim);
		framef_pi0_sim.setLocationRelativeTo(null);
		framef_pi0_sim.setVisible(true);
		canf_pi0_sim.divide(1, 1);
		canf_pi0_sim.cd(0);
		canf_pi0_sim.setFont("Arial");
		pi0_sim_hf.setTitle("f vs. #phi");
		pi0_sim_hf.setTitleX("#phi [deg]");
		pi0_sim_hf.setTitleY("f");
	//	pi0_sim_hf.setLineWidth(2);
		pi0_sim_hf.setMarkerSize(5);
		pi0_sim_hf.setMarkerStyle(1);
		canf_pi0_sim.getPad(0).setTitleFontSize(32);
		canf_pi0_sim.getPad(0).setAxisTitleFontSize(32);
		canf_pi0_sim.getPad(0).setAxisLabelFontSize(24);
		canf_pi0_sim.getPad(0).setAxisRange(0, 360, 0, 1.0);
		canf_pi0_sim.draw(pi0_sim_hf, "same");
		F1D ff_sim= new F1D("ff_sim", "[a]*(exp(-0.5((x-[m1])/[s])((x-[m1])/[s]))+exp(-0.5((x-[m2])/[s])((x-[m2])/[s]))) ", 0, 360);
		ff_sim.setParameter(0, 0.5);
		ff_sim.setParLimits(0, 0, 1);
		ff_sim.setParameter(1, 90);
		ff_sim.setParLimits(1, 65, 125);
		ff_sim.setParameter(2, 45);
		ff_sim.setParLimits(2, 10, 80);
		ff_sim.setParameter(3, 270);
		ff_sim.setParLimits(3, 235, 305);
		DataFitter.fit(ff_sim, pi0_sim_hf, "Q");
		ff_sim.setLineColor(2);
		ff_sim.setLineWidth(3);
		ff_sim.setOptStat(111110);
		canf_pi0_sim.draw(ff_sim, "same");
		double a = ff_sim.parameter(0).value();
		double da = ff_sim.parameter(0).error();
		double m1 = ff_sim.parameter(1).value();
		double dm1 = ff_sim.parameter(1).error();
		double s = ff_sim.parameter(2).value();
		double ds = ff_sim.parameter(2).error();
		double m2 = ff_sim.parameter(3).value();
		double dm2 = ff_sim.parameter(3).error();
		double e1 = Math.exp(-0.5*((90-m1)/s)*((90-m1)/s));
		double e2 = Math.exp(-0.5*((90-m2)/s)*((90-m2)/s));
		double f = a*(Math.exp(-0.5*((90-m1)/s)*((90-m1)/s))+Math.exp(-0.5*((90-m2)/s)*((90-m2)/s)));
		double ferr = Math.sqrt(((e1+e2)*da)*((e1+e2)*da)+(a*((90-m1)/s)*e1*dm1)*(a*((90-m1)/s)*e1*dm1)+(a*((90-m2)/s)*e2*dm2)*(a*((90-m2)/s)*e2*dm2)
							+(a*(((90-m1)/s)*e1+((90-m2)/s)*e2)*ds/s)*(a*(((90-m1)/s)*e1+((90-m2)/s)*e2)*ds/s));
//		pi0_sim_hf.addPoint(90, f, 0, ferr);
//		canf_pi0_sim.draw(pi0_sim_hf, "same");
//		System.out.println("A: " + fBSA_corr_pi0_sim.parameter(0).value() + " +/- " + fBSA_corr_pi0_sim.parameter(0).error());
//		System.out.println("B: " + fBSA_corr_pi0_sim.parameter(1).value() + " +/- " + fBSA_corr_pi0_sim.parameter(1).error());
		
		JFrame framebsa_corr_pi0_sim = new JFrame("#pi^0_sim Corrected BSA");
		framebsa_corr_pi0_sim.setSize(500, 500);
		EmbeddedCanvas canbsa_corr_pi0_sim = new EmbeddedCanvas();
		framebsa_corr_pi0_sim.add(canbsa_corr_pi0_sim);
		framebsa_corr_pi0_sim.setLocationRelativeTo(null);
		framebsa_corr_pi0_sim.setVisible(true);
		canbsa_corr_pi0_sim.cd(0);
		canbsa_corr_pi0_sim.setFont("Arial");
		corr_pi0_sim_hBSAy.setTitle("#pi^0_sim Corrected Beam Spin Asymmetry");
		corr_pi0_sim_hBSAy.setTitleX("#phi [deg]");
		corr_pi0_sim_hBSAy.setTitleY("#pi^0_sim Corrected BSA");
		corr_pi0_sim_hBSAy.setMarkerSize(5);
		corr_pi0_sim_hBSAy.setMarkerStyle(1);
		corr_pi0_sim_hBSAy.setLineThickness(1);
		canbsa_corr_pi0_sim.getPad(0).setTitleFontSize(32);
		canbsa_corr_pi0_sim.getPad(0).setAxisTitleFontSize(32);
		canbsa_corr_pi0_sim.getPad(0).setAxisLabelFontSize(24);
		canbsa_corr_pi0_sim.getPad(0).setAxisRange(0, 360, -0.3, 0.3);
		canbsa_corr_pi0_sim.getPad(0).setStatBoxFontSize(18);
		canbsa_corr_pi0_sim.draw(corr_pi0_sim_hBSAy, "same");
		F1D fBSA_corr_pi0_sim= new F1D("fBSA_corr_pi0_sim", "[A]*sin(x/57.3)/(1+[B]*cos(x/57.3))", 0, 360);
		fBSA_corr_pi0_sim.setParameter(0, 0.2);
		fBSA_corr_pi0_sim.setParameter(1, -0.5);
		DataFitter.fit(fBSA_corr_pi0_sim, corr_pi0_sim_hBSAy, "Q");
		fBSA_corr_pi0_sim.setLineColor(2);
		fBSA_corr_pi0_sim.setLineWidth(3);
		fBSA_corr_pi0_sim.setOptStat(1110);
		canbsa_corr_pi0_sim.draw(fBSA_corr_pi0_sim, "same");
		System.out.println("==========================================================================================");
		System.out.println("BSA_corr_pi0_sim:	x	#deltax	y	#delta	y");
		for(int p = 0; p < 30; p++)		
		{
			System.out.println("	" + phibc[p] + "	" + 0 + "	" + corr_pi0_sim_bsa.get(p) + "	" + corr_pi0_sim_bsaerr.get(p));
		}
		System.out.println("A: " + fBSA_corr_pi0_sim.parameter(0).value() + " +/- " + fBSA_corr_pi0_sim.parameter(0).error());
		System.out.println("B: " + fBSA_corr_pi0_sim.parameter(1).value() + " +/- " + fBSA_corr_pi0_sim.parameter(1).error());
		
		JFrame framer_pi0_sim_bin = new JFrame("Binned #pi^0_sim r");
		framer_pi0_sim_bin.setSize(1200, 600);
		EmbeddedCanvas canr_pi0_sim_bin = new EmbeddedCanvas();
		framer_pi0_sim_bin.add(canr_pi0_sim_bin);
		framer_pi0_sim_bin.setLocationRelativeTo(null);
		framer_pi0_sim_bin.setVisible(true);
		canr_pi0_sim_bin.divide(4, 2);
		for(int ih = 0; ih  < 8; ih++)
		{
			canr_pi0_sim_bin.cd(ih);
			canr_pi0_sim_bin.setFont("Arial");
			pi0_sim_hr_bin.getItem(ih).setMarkerSize(3);
			pi0_sim_hr_bin.getItem(ih).setMarkerStyle(1);
			pi0_sim_hr_bin.getItem(ih).setLineThickness(1);
			pi0_sim_hr_bin.getItem(ih).setTitle("Bin " + (ih+1));
			pi0_sim_hr_bin.getItem(ih).setTitleX("#phi [#degree]");
			pi0_sim_hr_bin.getItem(ih).setTitleY("#pi^0_sim r");
			canr_pi0_sim_bin.getPad(ih).setAxisLabelFontSize(10);
			canr_pi0_sim_bin.getPad(ih).setAxisTitleFontSize(10);
			canr_pi0_sim_bin.getPad(ih).setTitleFontSize(10);
			canr_pi0_sim_bin.getPad(ih).setAxisRange(0, 360, 0, 2);
			canr_pi0_sim_bin.draw(pi0_sim_hr_bin.getItem(ih), "same");
		}
		
		JFrame framef_pi0_sim_bin = new JFrame("Binned #pi^0_sim f");
		framef_pi0_sim_bin.setSize(1200, 600);
		EmbeddedCanvas canf_pi0_sim_bin = new EmbeddedCanvas();
		framef_pi0_sim_bin.add(canf_pi0_sim_bin);
		framef_pi0_sim_bin.setLocationRelativeTo(null);
		framef_pi0_sim_bin.setVisible(true);
		canf_pi0_sim_bin.divide(4, 2);
		for(int ih = 0; ih  < 8; ih++)
		{
			canf_pi0_sim_bin.cd(ih);
			canf_pi0_sim_bin.setFont("Arial");
			pi0_sim_hf_bin.getItem(ih).setMarkerSize(3);
			pi0_sim_hf_bin.getItem(ih).setMarkerStyle(1);
			pi0_sim_hf_bin.getItem(ih).setLineThickness(1);
			pi0_sim_hf_bin.getItem(ih).setTitle("Bin " + (ih+1));
			pi0_sim_hf_bin.getItem(ih).setTitleX("#phi [#degree]");
			pi0_sim_hf_bin.getItem(ih).setTitleY("#pi^0_sim f");
			canf_pi0_sim_bin.getPad(ih).setAxisLabelFontSize(10);
			canf_pi0_sim_bin.getPad(ih).setAxisTitleFontSize(10);
			canf_pi0_sim_bin.getPad(ih).setTitleFontSize(10);
			canf_pi0_sim_bin.getPad(ih).setAxisRange(0, 360, 0, 1.0);
			canf_pi0_sim_bin.draw(pi0_sim_hf_bin.getItem(ih), "same");
		}
		
		JFrame framebsa_corr_pi0_sim_bin = new JFrame("Binned Corrected #pi^0_sim BSA");
		framebsa_corr_pi0_sim_bin.setSize(1200, 600);
		EmbeddedCanvas canbsa_corr_pi0_sim_bin = new EmbeddedCanvas();
		framebsa_corr_pi0_sim_bin.add(canbsa_corr_pi0_sim_bin);
		framebsa_corr_pi0_sim_bin.setLocationRelativeTo(null);
		framebsa_corr_pi0_sim_bin.setVisible(true);
		canbsa_corr_pi0_sim_bin.divide(4, 2);
		for(int ih = 0; ih  < 8; ih++)
		{
			canbsa_corr_pi0_sim_bin.cd(ih);
			canbsa_corr_pi0_sim_bin.setFont("Arial");
			corr_pi0_sim_hBSAy_bin.getItem(ih).setMarkerSize(3);
			corr_pi0_sim_hBSAy_bin.getItem(ih).setMarkerStyle(1);
			corr_pi0_sim_hBSAy_bin.getItem(ih).setLineThickness(1);
			corr_pi0_sim_hBSAy_bin.getItem(ih).setTitle("Bin " + (ih+1));
			corr_pi0_sim_hBSAy_bin.getItem(ih).setTitleX("#phi [#degree]");
			corr_pi0_sim_hBSAy_bin.getItem(ih).setTitleY("#pi^0_sim Corrected BSA");
			canbsa_corr_pi0_sim_bin.getPad(ih).setAxisLabelFontSize(10);
			canbsa_corr_pi0_sim_bin.getPad(ih).setAxisTitleFontSize(10);
			canbsa_corr_pi0_sim_bin.getPad(ih).setTitleFontSize(10);
			canbsa_corr_pi0_sim_bin.getPad(ih).setAxisRange(0, 360, -0.3, 0.3);
			canbsa_corr_pi0_sim_bin.draw(corr_pi0_sim_hBSAy_bin.getItem(ih), "same");
			F1D fBSA_bin= new F1D("fBSA_bin", "[A]*sin(x/57.3)/(1+[B]*cos(x/57.3))", 0, 360);
			fBSA_bin.setParameter(0, 0.2);
			fBSA_bin.setParLimits(0, 0, 1.0);
			fBSA_bin.setParameter(1, -0.5);
			DataFitter.fit(fBSA_bin, corr_pi0_sim_hBSAy_bin.getItem(ih), "Q");
			if(fBSA_bin.parameter(0).value() > 0.199 && fBSA_bin.parameter(0).value() <0.201)
			{
				fBSA_bin.setParameter(0, 0.5);
				fBSA_bin.setParLimits(0, 0, 1.0);
				fBSA_bin.setParameter(1, -0.5);
				DataFitter.fit(fBSA_bin, corr_pi0_sim_hBSAy_bin.getItem(ih), "Q");
			}
			if(fBSA.parameter(0).value() > 0.499 && fBSA_bin.parameter(0).value() < 0.501)
			{
				fBSA_bin.setParameter(0, 0.0);
				fBSA_bin.setParLimits(0, 0, 1.0);
				fBSA_bin.setParameter(1, -0.5);
				DataFitter.fit(fBSA_bin, corr_pi0_sim_hBSAy_bin.getItem(ih), "Q");
			}
			fBSA_bin.setLineColor(2);
			fBSA_bin.setLineWidth(2);
			fBSA_bin.setOptStat(1110);
			canbsa_corr_pi0_sim_bin.draw(fBSA_bin, "same");
	//		System.out.println("Bin No.: " + ih);
	//		System.out.println("A: " + fBSA_bin.parameter(0).value() + " +/- " + fBSA_bin.parameter(0).error());
	//		System.out.println("B: " + fBSA_bin.parameter(1).value() + " +/- " + fBSA_bin.parameter(1).error());
		}
		
		JFrame frame_alphat_corr_pi0_sim = new JFrame("#pi^0_sim Corrected #alpha vs. -t");
		frame_alphat_corr_pi0_sim.setSize(1000, 1000);
		EmbeddedCanvas can_alphat_corr_pi0_sim = new EmbeddedCanvas();
		frame_alphat_corr_pi0_sim.add(can_alphat_corr_pi0_sim);
		frame_alphat_corr_pi0_sim.setLocationRelativeTo(null);
		frame_alphat_corr_pi0_sim.setVisible(true);
		can_alphat_corr_pi0_sim.divide(2, 2);
		can_alphat_corr_pi0_sim.cd(0);
		can_alphat_corr_pi0_sim.setFont("Arial");
		corr_pi0_sim_halpha_vs_t.getItem(0).setTitle("#pi^0_sim Corrected BSA #alpha vs. -t");
		corr_pi0_sim_halpha_vs_t.getItem(0).setTitleX("-t [GeV^2]");
		corr_pi0_sim_halpha_vs_t.getItem(0).setTitleY("#pi^0_sim Corrected #alpha");
		corr_pi0_sim_halpha_vs_t.getItem(0).setMarkerSize(5);
		corr_pi0_sim_halpha_vs_t.getItem(0).setMarkerStyle(1);
		corr_pi0_sim_halpha_vs_t.getItem(0).setLineThickness(1);
		can_alphat_corr_pi0_sim.getPad(0).setTitleFontSize(32);
		can_alphat_corr_pi0_sim.getPad(0).setAxisTitleFontSize(32);
		can_alphat_corr_pi0_sim.getPad(0).setAxisLabelFontSize(24);
		can_alphat_corr_pi0_sim.getPad(0).setAxisRange(0, 2, -0.1, 0.45);
		can_alphat_corr_pi0_sim.getPad(0).setStatBoxFontSize(18);
		can_alphat_corr_pi0_sim.draw(corr_pi0_sim_halpha_vs_t.getItem(0), "same");
		can_alphat_corr_pi0_sim.cd(1);
		can_alphat_corr_pi0_sim.setFont("Arial");
		corr_pi0_sim_halpha_vs_t.getItem(1).setTitle("#pi^0_sim Corrected BSA #alpha vs. -t");
		corr_pi0_sim_halpha_vs_t.getItem(1).setTitleX("-t [GeV^2]");
		corr_pi0_sim_halpha_vs_t.getItem(1).setTitleY("#pi^0_sim Corrected #alpha");
		corr_pi0_sim_halpha_vs_t.getItem(1).setMarkerSize(5);
		corr_pi0_sim_halpha_vs_t.getItem(1).setMarkerStyle(1);
		corr_pi0_sim_halpha_vs_t.getItem(1).setLineThickness(1);
		can_alphat_corr_pi0_sim.getPad(1).setTitleFontSize(32);
		can_alphat_corr_pi0_sim.getPad(1).setAxisTitleFontSize(32);
		can_alphat_corr_pi0_sim.getPad(1).setAxisLabelFontSize(24);
		can_alphat_corr_pi0_sim.getPad(1).setAxisRange(0, 2, -0.1, 0.45);
		can_alphat_corr_pi0_sim.getPad(1).setStatBoxFontSize(18);
		can_alphat_corr_pi0_sim.draw(corr_pi0_sim_halpha_vs_t.getItem(1), "same");
		corr_pi0_sim_halpha_vs_t.getItem(2).setMarkerSize(5);
		corr_pi0_sim_halpha_vs_t.getItem(2).setLineThickness(1);
		corr_pi0_sim_halpha_vs_t.getItem(2).setMarkerStyle(0);
		corr_pi0_sim_halpha_vs_t.getItem(2).setMarkerColor(2);
		corr_pi0_sim_halpha_vs_t.getItem(2).setLineColor(2);
		can_alphat_corr_pi0_sim.draw(corr_pi0_sim_halpha_vs_t.getItem(2), "same");
		can_alphat_corr_pi0_sim.cd(2);
		can_alphat_corr_pi0_sim.setFont("Arial");
		corr_pi0_sim_halpha_vs_t.getItem(3).setTitle("#pi^0_sim Corrected BSA #alpha vs. -t");
		corr_pi0_sim_halpha_vs_t.getItem(3).setTitleX("-t [GeV^2]");
		corr_pi0_sim_halpha_vs_t.getItem(3).setTitleY("#pi^0_sim Corrected #alpha");
		corr_pi0_sim_halpha_vs_t.getItem(3).setMarkerSize(5);
		corr_pi0_sim_halpha_vs_t.getItem(3).setMarkerStyle(1);
		corr_pi0_sim_halpha_vs_t.getItem(3).setLineThickness(1);
		can_alphat_corr_pi0_sim.getPad(2).setTitleFontSize(32);
		can_alphat_corr_pi0_sim.getPad(2).setAxisTitleFontSize(32);
		can_alphat_corr_pi0_sim.getPad(2).setAxisLabelFontSize(24);
		can_alphat_corr_pi0_sim.getPad(2).setAxisRange(0, 2, -0.1, 0.45);
		can_alphat_corr_pi0_sim.getPad(2).setStatBoxFontSize(18);
		can_alphat_corr_pi0_sim.draw(corr_pi0_sim_halpha_vs_t.getItem(3), "same");
		corr_pi0_sim_halpha_vs_t.getItem(4).setMarkerSize(5);
		corr_pi0_sim_halpha_vs_t.getItem(4).setLineThickness(1);
		corr_pi0_sim_halpha_vs_t.getItem(4).setMarkerStyle(0);
		corr_pi0_sim_halpha_vs_t.getItem(4).setMarkerColor(2);
		corr_pi0_sim_halpha_vs_t.getItem(4).setLineColor(2);
		can_alphat_corr_pi0_sim.draw(corr_pi0_sim_halpha_vs_t.getItem(4), "same");
		corr_pi0_sim_halpha_vs_t.getItem(5).setMarkerSize(5);
		corr_pi0_sim_halpha_vs_t.getItem(5).setLineThickness(1);
		corr_pi0_sim_halpha_vs_t.getItem(5).setMarkerStyle(2);
		corr_pi0_sim_halpha_vs_t.getItem(5).setMarkerColor(4);
		corr_pi0_sim_halpha_vs_t.getItem(5).setLineColor(4);
		can_alphat_corr_pi0_sim.draw(corr_pi0_sim_halpha_vs_t.getItem(5), "same");
		can_alphat_corr_pi0_sim.cd(3);
		can_alphat_corr_pi0_sim.setFont("Arial");
		corr_pi0_sim_halpha_vs_t.getItem(6).setTitle("#pi^0_sim Corrected BSA #alpha vs. -t");
		corr_pi0_sim_halpha_vs_t.getItem(6).setTitleX("-t [GeV^2]");
		corr_pi0_sim_halpha_vs_t.getItem(6).setTitleY("#pi^0_sim Corrected #alpha");
		corr_pi0_sim_halpha_vs_t.getItem(6).setMarkerSize(5);
		corr_pi0_sim_halpha_vs_t.getItem(6).setMarkerStyle(1);
		corr_pi0_sim_halpha_vs_t.getItem(6).setLineThickness(1);
		can_alphat_corr_pi0_sim.getPad(3).setTitleFontSize(32);
		can_alphat_corr_pi0_sim.getPad(3).setAxisTitleFontSize(32);
		can_alphat_corr_pi0_sim.getPad(3).setAxisLabelFontSize(24);
		can_alphat_corr_pi0_sim.getPad(3).setAxisRange(0, 2, -0.1, 0.45);
		can_alphat_corr_pi0_sim.getPad(3).setStatBoxFontSize(18);
		can_alphat_corr_pi0_sim.draw(corr_pi0_sim_halpha_vs_t.getItem(6), "same");
		corr_pi0_sim_halpha_vs_t.getItem(7).setMarkerSize(5);
		corr_pi0_sim_halpha_vs_t.getItem(7).setLineThickness(1);
		corr_pi0_sim_halpha_vs_t.getItem(7).setMarkerStyle(0);
		corr_pi0_sim_halpha_vs_t.getItem(7).setMarkerColor(2);
		corr_pi0_sim_halpha_vs_t.getItem(7).setLineColor(2);
		can_alphat_corr_pi0_sim.draw(corr_pi0_sim_halpha_vs_t.getItem(7), "same");
		
	}
}