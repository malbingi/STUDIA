// moments
        double m_rho, m_e, m_eps, m_jx, m_jy, m_qx, m_qy, m_pxx, m_pxy;
        double m_rhoe, m_ee, m_epse, m_jxe, m_jye, m_qxe, m_qye, m_pxxe, m_pxye;
        m_rho = rho;
        m_jx = rho * ux;
        m_jy = rho * uy;
        m_e = -(4.0 * fc + fe + fn + fw + fs) + 2.0 * (fne + fnw + fse + fsw);
        m_eps = (4.0 * fc + fne + fnw + fse + fsw) - 2.0 * (fe + fn + fw + fs);
        m_qx = 2.0 * (fw - fe) + fne - fnw - fsw + fse;
        m_qy = 2.0 * (fs - fn) + fne + fnw - fsw - fse;
        m_pxx = fe - fn + fw - fs;
        m_pxy = fne - fnw + fsw - fse;
      // equilibria
        m_rhoe = rho;
        m_jxe = rho * ux;
        m_jye = rho * uy;

        m_ee = rho * (-2.0 + 3.0 * (ux * ux + uy * uy));
        m_epse = rho * (1.0 - 3.0 * (ux * ux + uy * uy));
        m_qxe = -rho * ux;
        m_qye = -rho * uy;
        m_pxxe = rho * (ux * ux - uy * uy);
        m_pxye = rho * ux * uy;

        double om_e, om_eps, om_q, om_nu;
        om_e = 1.63;
        om_eps = 1.54;
        om_q = 1.0;
        om_nu = 1.0 / tau;

        //collision

        m_e += om_e * (m_ee - m_e);
        m_eps += om_eps * (m_epse - m_eps);
        m_qx += om_q * (m_qxe - m_qx);
        m_qy += om_q * (m_qye - m_qy);
        m_pxx += om_nu * (m_pxxe - m_pxx);
        m_pxy += om_nu * (m_pxye - m_pxy);

        fstar[_C] = 1. / 9. * (m_rho - m_e + m_eps);
        fstar[_E] = 1. / 36. * (4. * m_rho - m_e - 2. * m_eps) + 1. / 6. * (m_jx - m_qx) + m_pxx * .25;
        fstar[_N] = 1. / 36. * (4. * m_rho - m_e - 2. * m_eps) + 1. / 6. * (m_jy - m_qy) - m_pxx * .25;
        fstar[_W] = 1. / 36. * (4. * m_rho - m_e - 2. * m_eps) + 1. / 6. * (-m_jx + m_qx) + m_pxx * .25;
        fstar[_S] = 1. / 36. * (4. * m_rho - m_e - 2. * m_eps) + 1. / 6. * (-m_jy + m_qy) - m_pxx * .25;
        fstar[_NE] = 1. / 36. * (4. * m_rho + 2. * m_e + m_eps) + 1. / 12. * (2. * m_jx + 2. * m_jy + (m_qx + m_qy)) + m_pxy * .25;   
        fstar[_NW] = 1. / 36. * (4. * m_rho + 2. * m_e + m_eps) + 1. / 12. * (-2. * m_jx + 2. * m_jy + (-m_qx + m_qy)) - m_pxy * .25; 
        fstar[_SW] = 1. / 36. * (4. * m_rho + 2. * m_e + m_eps) - 1. / 12. * (2. * m_jx + 2. * m_jy + (m_qx + m_qy)) + m_pxy * .25;   
        fstar[_SE] = 1. / 36. * (4. * m_rho + 2. * m_e + m_eps) + 1. / 12. * (2. * m_jx - 2. * m_jy + (m_qx - m_qy)) - m_pxy * .25;   
