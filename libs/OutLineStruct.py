#!/usr/bin/env python

MULTIPLEGENOTYPE = "multiple_genotype"
NOTENOUGHVARIANTS = "var_number"

class OutLineStruct:
    def __init__(self, name, pos, ref, var, ref_num, var_num,lh, 
                so, var_all, bias, gq, pl, gt_num, gt_class):
        # type: (str, int, str, str, int, int, float, float, float, float, float)
        self.name = name
        self.pos = pos
        self.ref = ref
        self.var = var
        self.ref_num = ref_num
        self.var_num = var_num
        self.total_num = ref_num + var_num
        self.ad = '{},{}'.format(ref_num, var_num)
        self.lh = lh # [Seq.  Noise, Ampl. Artefact/Error, Het. SNV, Hom. SNV]
        self.so = so
        self.var_all = var_all
        self.bias = bias
        self.gq = gq
        self.pl = pl
        self.get = ''
        self.gt_num = gt_num
        self.gt_class = gt_class


    def __str__(self):
        return '{}\t{}\t{}\t{}\t{}\t{}\t{}' \
            .format(self.name, self.pos, self.ref, self.var, self.ref_num,
                self.var_num, self.so)


    def get_varcall_str(self):
        return '{0}\t{1}\t{1}\t{2}\t{3}\t{4}\t{5}\tPASS\t{6}\t{7}\t{8}\t{9}\t{10}' \
            .format(self.name, self.pos, self.ref, self.var, self.ref_num,
                self.var_num, self.bias, self.lh[0], self.lh[1], self.lh[2],
                self.lh[3])


    def get_reason_str(self):
        return '{}\t{}\t{}\t{}\t{}\t{}\t{}' \
            .format(self.name, self.pos, self.ref, self.var, self.ref_num,
                self.var_num, self.so)


    def get_vcf_str(self, eta, min_var, min_var_frac):
        gt = self._get_gt(eta, min_var_frac)
        if self.so != "":
            format_str = "GT:SO:BN:AD:BI:GQ:FPL"
            cell_str = "{}:{}:{}:{}:{:.3f}:{}:{}" \
                .format(gt, self._transform_so(), self.so, self.ad, self.bias,
                    self.gq, self.pl)
        else:
            format_str = "GT:AD:BI:GQ:FPL"
            cell_str = "{}:{}:{:.3f}:{}:{}" \
                .format(gt, self.ad, self.bias, self.gq, self.pl)

        filter_str = self._get_filter(min_var)
        
        return '{}\t{}\t.\t{}\t{}\t{}\t{}\tNS=1\t{}\t{}' \
            .format(self.name, self.pos, self.ref, self.var_all,
                self.gq, filter_str, format_str, cell_str)


    def _get_gt(self, eta, min_var_frac):
        # <TODO, NB> Why 2 * Sequencing  Noise and Ampl. Artefact/Error / eta
        lls = [2 * self.lh[0], self.lh[1] / eta, self.lh[2], self.lh[3]]
        max_ll_id = lls.index(max(lls))
        if max_ll_id < 2:
            result_str = "0/0"
        elif max_ll_id == 2:
            # <DONE, NB>
            if self.var_num >= (self.var_num + self.ref_num) * min_var_frac:
            # <BUG, Original> Why: 7 * var > ref ???
            # if 7 * self.var_num > self.ref_num:
                result_str = "0/1"
            else:
                result_str = "0/0"
        else:
            result_str = "1/1"

        if self.gt_class == 2:
            if result_str == "0/0":
                result_str = "1/1"
            elif result_str == "0/1":
                result_str = "1/2"
            else:
                result_str = "1/1"

        return result_str


    def _get_filter(self, min_var):    
        if self.gt_class == 1:
            var_num = self.var_num
        else:
            var_num = self.ref_num

        if len(self.var.split(",")) < self.gt_num or var_num < min_var:
            if var_num < min_var:
                comment = '{},{}'.format(MULTIPLEGENOTYPE, NOTENOUGHVARIANTS)
            else:
                comment = MULTIPLEGENOTYPE
        else:
            comment = "."
        return comment


    def _transform_so(self):
        if self.so == "ref":
            return "True"
        elif self.so.startswith("var"):
            return "False"
        else:
            return "NA"


if __name__ == '__main__':
    print('Here be dragons...')