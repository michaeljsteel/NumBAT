
#ifndef _numbat_decl_h_
#define _numbat_decl_h_

#define RETONERROR(ec) if (ec .ne. 0) then ; return ; endif

#define RET_ON_NBERR(nberr) if (nberr%errco .ne. 0) then ; return ; endif
#define RET_ON_NBERR_UNFOLD(nberr) if (nberr%errco .ne. 0) then ; call nberr%to_py(errco, emsg); return ; endif

#define NBERROR_101 101_8
#define NBERROR_102 102_8
#define NBERROR_103 103_8
#define NBERROR_104 104_8
#define NBERROR_105 105_8
#define NBERROR_106 106_8
#define NBERROR_107 107_8
#define NBERROR_108 108_8
#define NBERROR_109 109_8
#define NBERROR_110 110_8
#define NBERROR_111 111_8
#define NBERROR_112 112_8
#define NBERROR_113 113_8
#define NBERROR_114 114_8
#define NBERROR_115 115_8
#define NBERROR_116 116_8
#define NBERROR_117 117_8
#define NBERROR_118 118_8
#define NBERROR_119 119_8
#define NBERROR_120 120_8
#define NBERROR_121 121_8
#define NBERROR_122 122_8
#define NBERROR_123 123_8
#define NBERROR_124 124_8
#define NBERROR_125 125_8
#define NBERROR_126 126_8
#define NBERROR_127 127_8
#define NBERROR_128 128_8
#define NBERROR_129 129_8
#define NBERROR_130 130_8






#endif

