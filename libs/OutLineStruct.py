#!/usr/bin/env python

import numpy as np

MULTIPLEGENOTYPE = 'multiple_genotype'
NOTENOUGHVARIANTS = 'var_number'

class OutLineStruct:
    def __init__(self, name, pos, ref, var, ref_num, var_num, log_lh, 
                so, var_all, bias, gt_num, gt_class):
        # type: (str, int, str, str, int, int, float, float, float, float, float)
        self.name = name
        self.pos = pos
        self.ref = ref
        self.var = var
        self.ref_num = ref_num
        self.var_num = var_num
        self.total_num = ref_num + var_num
        self.ad = '{},{}'.format(ref_num, var_num)
        self.log_lh = log_lh # [REF/REF, REF/ALT, ALT/ALT]
        self.so = so
        self.var_all = var_all
        self.bias = bias
        self.gt_num = gt_num
        self.gt_class = gt_class


    def __str__(self):
        return '{}\t{}\t{}\t{}\t{}\t{}\t{}' \
            .format(self.name, self.pos, self.ref, self.var, self.ref_num,
                self.var_num, self.so)


    def get_varcall_str(self):
        return '{0}\t{1}\t{1}\t{2}\t{3}\t{4}\t{5}\tPASS\t{6}\t{7}\t{8}\t{9}' \
            .format(self.name, self.pos, self.ref, self.var, self.ref_num,
                self.var_num, self.bias, self.log_lh[0], self.log_lh[1], 
                self.log_lh[2])


    def get_reason_str(self):
        return '{}\t{}\t{}\t{}\t{}\t{}\t{}' \
            .format(self.name, self.pos, self.ref, self.var, self.ref_num,
                self.var_num, self.so)


    def get_vcf_str(self, eta, min_var, min_var_frac):
        gt = self._get_gt(eta, min_var_frac)
        pl = self._get_pl()
        gq = self._get_gq(pl)
       
        if self.so != '':
            format_str = 'GT:SO:BN:AD:BI:GQ:PL'
            cell_str = '{}:{}:{}:{}:{:.3f}:{}:{:.0f},{:.0f},{:.0f}' \
                .format(gt, self._transform_so(), self.so, self.ad, self.bias,
                    gq, *pl)
        else:
            format_str = 'GT:AD:BI:GQ:PL'
            cell_str = '{}:{}:{:.3f}:{}:{:.0f},{:.0f},{:.0f}' \
                .format(gt, self.ad, self.bias, gq, *pl)

        filter_str = self._get_filter(min_var)
        
        return '{}\t{}\t.\t{}\t{}\t{}\t{}\tNS=1\t{}\t{}' \
            .format(self.name, self.pos, self.ref, self.var_all, gq, 
                filter_str, format_str, cell_str)


    def _get_gt(self, eta, min_var_frac):
        # <TODO, NB> Why 2 * Sequencing  Noise and Ampl. Artefact/Error / eta
        # lls = [2 * self.lh[0], self.lh[1] / eta, self.lh[2], self.lh[3]]
        max_ll_id = np.argmax(self.log_lh)
        if max_ll_id == 0:
            result_str = '0/0'
        elif max_ll_id == 1:
            # <DONE, NB>
            # Apply 'Likelihood ratio test' only to test for H_0/H_1
            if self.var_num >= (self.var_num + self.ref_num) * min_var_frac \
                    and self.log_lh[1] >= self.log_lh[0] - np.log10(eta):
            # <BUG, Original> Why: 7 * var > ref ???
            # if 7 * self.var_num > self.ref_num:
                result_str = '0/1'
            else:
                result_str = '0/0'
        else:
            result_str = '1/1'

        if self.gt_class == 2:
            if result_str == '0/0':
                result_str = '1/1'
            elif result_str == '0/1':
                result_str = '1/2'
            else:
                result_str = '1/1'

        return result_str


    def _get_pl(self):
        pl_raw = -10 * self.log_lh
        # <DONE, NB>
        pl = pl_raw - pl_raw.min()
        # <BUG, Original> Dont normalize Phred scores to smallest = 0
        # pl = pl_raw
        return pl


    # <INFO, NB> Calculation of GQ value (phred-scaled genotype quality score)
    # Replaced with np array for clarity and performance
    def _get_gq(self, pl):
        # type: (float, float, float) -> str
        value = sorted(pl)[1]
        return int(min(99, round(value)))


    def _get_filter(self, min_var):    
        if self.gt_class == 1:
            var_num = self.var_num
        else:
            var_num = self.ref_num

        if len(self.var.split(',')) < self.gt_num or var_num < min_var:
            if var_num < min_var:
                comment = '{};{}'.format(MULTIPLEGENOTYPE, NOTENOUGHVARIANTS)
            else:
                comment = MULTIPLEGENOTYPE
        else:
            comment = '.'
        return comment


    def _transform_so(self):
        if self.so == 'ref':
            return 'True'
        elif self.so.startswith('var'):
            return 'False'
        else:
            return 'NA'


if __name__ == '__main__':
    print('Here be dragons...')