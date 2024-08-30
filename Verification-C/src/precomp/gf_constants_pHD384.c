#include <stddef.h>
#include <stdint.h>
#include <tutil.h>
#include <fp2.h>
#include "gf_constants.h"

const fp2_t NQR_TABLE[20] = { {{ 0xcdf84757f0b5566, 0xb065e9dce582064d, 0x374193b0a8d289f2, 0x083c1e98b3d19a32, 0xdc269c89a54888e2, 0x1a28cf4a0cd7dec7}, {0x43d473ffd01bf4, 0xe7c6815e55bb04ae, 0xd270ad90d1c3c16c, 0xf1136bb9c7162794, 0xbae75ab3855e3237, 0x55748ac4420b494d}}, {{ 0x309887b0edee3a, 0x8d407960de403c50, 0x43fd415f083d0f18, 0xbd0bf2c1d1629a74, 0xfad2fa7085a313d8, 0x2ca7123e9ccccaae}, {0x10c447e922fbe1d, 0x7495fc6c31bd1894, 0xa4241cb06a7acff2, 0xa8a0f907e5a85176, 0x4a44870f4d4c09c5, 0x16345cca2c42a76d}}, {{ 0x2eba9dc6fa2deac, 0x687eeef79d385ff7, 0xd60cb27eaf5ff6e5, 0x56085b5c84cc1e34, 0x7e9fd4a973ec6983, 0x21978589b1bde625}, {0xd28b88079211b9, 0xccc4bf3ebaa49320, 0x723b15a1a8011fcd, 0xca811f8e857f7eae, 0x7acd42f1ded989b3, 0x80b04006cb112090}}, {{ 0x4532e84d02051d, 0x83080a25a87fe5cc, 0x1cedc4426121db4e, 0x0a7d66dcdd7d7006, 0x58e948ef53a19f5d, 0x2b4c41807883da93}, {0xd5afc4f83cede8, 0x1a55e30b85766593, 0x812a6e9c1f9d0089, 0x6fcedc6356bd34fe, 0x32cbf501ecde462a, 0xfc0f1e01b3d5060c}}, {{ 0x28e5a047c413900, 0x4be1fb2999886151, 0x1489f3350388e7d2, 0xa67f90da6e307994, 0x8704ca68dbd1103f, 0x2dedaf5478fbcb47}, {0xc39f0abc58d70eb, 0x1351d5a233b037b4, 0x6fe89ebf217ea655, 0xf4b25f10c4b29fad, 0x72bd163955f4befa, 0x1c7ac905a20a8c65}}, {{ 0x3ec0092754ea28e, 0x3838d5cf27c771bf, 0x0541bc008f7c8973, 0xabbdbecdc5bbe9a6, 0x1da8fdf8e6941489, 0x1b89fd3ac5c1bedc}, {0xc1bcc249d1e79ee, 0x2ba640462e7c2408, 0x482a5b419fcbc307, 0xd1c995237cb44099, 0xb9fbad59fc39d224, 0x1dce73cf7f59ff33}}, {{ 0x1db101650d582fb, 0xf1eb82321bde30d8, 0xef6a89ecfdff52ca, 0x355cc7374c1ca42c, 0x9903dfc94e480b8b, 0x33a4af74ac8072c4}, {0xa11663d0d6b8ebf, 0x16dc1d7b8512656d, 0x8c273f52706560aa, 0x22c5cbd26ff6a2e8, 0xf4cb00caa909b911, 0x3ff8dd0571b9f912}}, {{ 0xa35f46a76e4702c, 0xf72b3057022352d4, 0x812aa5d8f5a5773c, 0xba4c1e978ddb949a, 0x5c88c70a170b11c5, 0x13c26b9e9728b0fe}, {0xea41da8f5010381, 0x772122246d3d99f9, 0x15df98dc00da9667, 0xa6d03f9ca3051765, 0xbd4f40a7d80e3fc5, 0x1cb443d28c2be3fd}}, {{ 0x238d076537a4301, 0xd2cb9879afadf8ed, 0x21861d996871cb0b, 0x6ea3315aa9cfa06c, 0x48a913eb1eca1460, 0x40b0fd520a4f1798}, {0xf12451c79fac0ed, 0x161b1203f925d693, 0xb8eeca4063f42320, 0x8e3c13f2aa8999fc, 0xfeb1336027caff5f, 0x2198c364b58bf5c1}}, {{ 0xc29226f423a3648, 0x918651f3f137b2fa, 0x4d9691c9b3998af4, 0x906c806519afd9fc, 0x4c93ade22baf506c, 0x3ea9ed2ff6b3547e}, {0x76a81bbf3ec3f89, 0x87d8dfb15551491a, 0x9d953a73b1330bbf, 0x6228c490623a7cd7, 0x218b241514132c6d, 0x24150d74fc77ee36}}, {{ 0x0c65ea781502ca4, 0x661c480459c876cd, 0x8fa97d23e83b2afd, 0xcd45b98655bddd62, 0xead80e1534065956, 0x15567746f2de68a6}, {0x50c7cbe862fe54b, 0x606fa28edf77e53d, 0x949390ca8602f8fc, 0x3ef68870e3a6cc04, 0xc4bd16901fe438f0, 0x1145ec0ef4c52d04}}, {{ 0xa82997172700acd, 0x3527e2163b75c41d, 0xf08aa77fbd7df789, 0xb888ef13fd7724f7, 0x0be004334dea2a39, 0x1a239f266831d4df}, {0xb3bf433a48b8629, 0xd6a506207024f0f4, 0xa2d41a9635cc956f, 0x27607a0c604a4104, 0xeeebd042a50cb1c7, 0x3dba637b55b3a4c3}}, {{ 0xd355b5575c4c600, 0xe7acd8d760cff5d8, 0x9dca514ce1f64d6b, 0xb3b3342e066a72e0, 0x7e4a5c03908a2430, 0x2e9506d88a5a68d0}, {0x182eac11bd5c5f9, 0xe84b7c219accf9eb, 0x30cd3abc35c469bb, 0x7c9da2991d9d8c64, 0xf6c3a7225f87f2b0, 0x105846581d428d31}}, {{ 0xa115052f48f7bc, 0x3f96b617ffd0ce0f, 0x225c1ead377f9055, 0x517b6f31836eb1e3, 0xca552f64ee2e6a4c, 0xc6b1b508bbd0b425}, {0x6154acb42dd964f, 0x048a4742e08649b0, 0xed81b80abeb72746, 0x480067a76c5fc121, 0xe88e082c2bc24852, 0x232c9ee53ec12e73}}, {{ 0x8a5e142cb56b02, 0x2ae5d5404effb6b8, 0x23f0acf79ea856f4, 0x3a988be641a88d58, 0x6447fd2f4c28e80b, 0x86e339678792a4fb}, {0xec5b7c6386cd77f, 0x6b6185c07ed98a01, 0x0c4deb7eb6ed23a4, 0x76be6e1ff2136404, 0x755db4195826a7e3, 0x1b36582f22502c6d}}, {{ 0xd1b3dab2ab23df8, 0xad1be25acbc32a7d, 0x53aa59b101dc0edb, 0xd778b6396106e890, 0x770fb5b9263d393c, 0x289f4081f9318161}, {0xc87e50e5d10ad87, 0xf6775711cf310250, 0xb55f6ac540cb96ec, 0xd3354c36dbf16cd6, 0x3d5cdf66432a5c1f, 0x399136b5791dde8c}}, {{ 0x6c1947575b05b4a, 0x310f9f8d370661e3, 0x5acdb49a8b464245, 0x9aff87a81b620520, 0x7799dfda810e321f, 0x11e1f4aafd7f95ab}, {0x76222af28f1c3b8, 0xb2b60abb9b624f10, 0xbd4b289f6cf436ea, 0x515a58bcd6ff7c9e, 0x433b288357a122b1, 0x114dbe594de260b5}}, {{ 0xf4c40eb2a297aff, 0xd3bf8ffbbe5d7ffe, 0xde2f626171f405f8, 0x852d9e52e5c30ab6, 0x639def6e9e148537, 0x38a8ad3d1686c35b}, {0x7973cd343c7c66b, 0x969bc7d062480a90, 0x415ee5aa54f06151, 0xfc70feb07aa38c26, 0x5150b1e1b244d3e0, 0x276ad78486633215}}, {{ 0x0075c0841110b66, 0x38970e29c189ce9a, 0x2958539486d0478c, 0x83a212279933c607, 0x04bf15150270fd0f, 0x22cca0f17207b4c0}, {0x7928bcdc1b0451, 0x71ad4c2d7ac149a4, 0x4c5dc477ff8bb3e7, 0x336b6dafd7ab7ec8, 0xd3d0cab4f447cf2e, 0x9c6a4dab757e048e}}, {{ 0x3cbacd25a399e72, 0x60acd591c97bfe88, 0x81b9b43bb86c62a7, 0x69981d2183957309, 0xaac673401aa91242, 0x2e14734432204104}, {0x55d3ca9918cf67, 0x603fb9f0a922d609, 0x3ee1276da5a80fc6, 0x09b05d1d15a3a7a3, 0xa87a6c77421452a6, 0x28fe8d1851e06a6b}}};
const fp2_t Z_NQR_TABLE[20] = { {{ 0x0c4d452a15cac0, 0x60781a092fec9257, 0x60f5633bfa00cb10, 0x69067c950eea9b5a, 0xaebda7387ee595e1, 0x173feedc06f36a7b}, {0xfa011cc08d314ae, 0x51f7accaed34c61c, 0xe10f328eea7b6403, 0x6bbfc86089574d06, 0xd1915d5305733e06, 0x31ec4e9b6ae3588d}}, {{ 0xf63f2431db218e, 0x1174c26170150cfb, 0x146f74af31ef926c, 0x3ddff7894f5aefb6, 0xbd140fc13367e469, 0xa165fc1e2be3a25e}, {0xc749ae3b0d2019, 0x04942fa8c41a2483, 0xdfca92daf86dbca5, 0xc2fb9148c8ad3727, 0x9e670f9c68b74d93, 0xafee116f334ba54a}}, {{ 0xe54e25496d86d0b, 0x0c548477c29591f7, 0xf91ff5a4f50522bc, 0xe025c3120b7546e0, 0x0502af8a47340875, 0x1d8d0815bb28943b}, {0x9d82ae1cd4be56, 0xb68ae2030341f3a1, 0x85b12771f68a1c49, 0xb7fe12ab8abc37b8, 0xa1ddbf16eb5d1cc8, 0xfe4216d1aaeeddd0}}, {{ 0xcaf8bfc515ccdce, 0xb57837d2b32d991b, 0xd02c3c1dd985a9fd, 0xa016eeff5c172eb3, 0xe0579173af872823, 0x388da98c4a3b0b7d}, {0xbb7c3080a69ca6e, 0x11466c81775d64f9, 0xb75425d5aab92bb2, 0x81859c7f65108445, 0x5ae3fbb6c9385ff1, 0x2b3f8b9995fe06f6}}, {{ 0x1b174551c492e25, 0x67fc7fe8361c1d0e, 0x2964c0d357b640d1, 0x84f1d9cd7c5d7970, 0x3be448cb09565cd9, 0x1bf250dc960dfebf}, {0x59f5a4e4973b4a, 0x23514b0aaa44eb64, 0x7bcdaf889995c05e, 0x5870c91cdc6697aa, 0xe428bd350500621f, 0x40627f339bbaf8bf}}, {{ 0x9b058f3e5453ef9, 0x8878a7a53b04d5a1, 0x739c8c02ccdc62c2, 0x94ef0794c8a25358, 0xe7c40924cc99cdd2, 0x20c79c4a92522923}, {0x30400575dfd154, 0xf06dd2bde70ce93a, 0x24b0281cefa47c12, 0xf0fdf1b3cdaf0972, 0x585148b9c460f927, 0xdf11b56a50b53513}}, {{ 0xfa3be09b4ec73c3, 0xce189806b7eecc38, 0x0411b538665ac3ec, 0x15c5e674a3e783ec, 0xf6a3533df8f42c0c, 0x21a9ba58030b2f5c}, {0x217ff0ba7271df, 0x744c6ac9021d6f2e, 0x87731c8329b3fe7e, 0x39148596a69f5ed9, 0x848e1acc479b11ee, 0xd8a2bea7c9a01cc5}}, {{ 0x57566635080ad1, 0x24032d99ce90e626, 0x8b0190d1b7ca303e, 0x95fb2db5bf2a58c1, 0xdcf75d1a1777fe0a, 0xffd42ac23177eaf6}, {0x3f0843fa9c17daf, 0xf67649a319ae7c28, 0xd8026679b8b3426c, 0x6cdd223f2375cb13, 0x22d504ab32da8350, 0x3617596b9ac0142a}}, {{ 0x815a8b6005a81e1, 0xdb0b491b54f9b703, 0x319af179a7fb29ae, 0x3083ffb6de94bb76, 0x46e9de366ae276d6, 0x17c9413c1fbe392a}, {0x3003da15f946d0f, 0xc7cb28de37697306, 0x8ef671d574bd4a0f, 0x9d1322d8ca74d671, 0x96dc3e47707b4bb6, 0x314b9cac6185fd6b}}, {{ 0xdcc72e4239d6949, 0x7a45bcca3de81f63, 0xf9757d80b7382b18, 0xa48032a4b8b9e8dc, 0x917b34f7cacb4ee0, 0x2391bafdb2482e2d}, {0x45053f80706f515, 0x1bc81b8a5345c76f, 0x41d2bd41062dbaee, 0x70149b470e848e82, 0x36655da7ea3098a8, 0x341f062871815648}}, {{ 0xb0513a487d3e42c, 0xa6512c517a93c92f, 0x78da639f91d2d4d3, 0x3e491948c3f89a7d, 0x5d99d7962f001f5e, 0x201b756839cd46f4}, {0xd320b343f09ab, 0x6846dfbc8d2656dc, 0x43f779dc94f041ab, 0xa75c178428707788, 0xd76eacce2bd14fd1, 0x463fc66ee53492b5}}, {{ 0x1140f2f84a4bed6, 0xa156b6568c2eb387, 0x9006877871a76877, 0x2151d2933d650705, 0xef632229f6e8d19a, 0x2d8bc3ae0ff28ff5}, {0x244c53b9372729, 0xa1c9ee2b5dcc1ece, 0xf5ea603a2e8a5ac8, 0xa548f746bbcc06f8, 0x06726f6a1fe4d3fa, 0xb5789694d66b156d}}, {{ 0x4f881210fe5e506, 0x5ca507c8749f73c1, 0x681d59cfd5e30a85, 0xaaf40898800a8489, 0xf71dac59e04f46d5, 0x34822f63715e1f09}, {0xa009f172c78e19e, 0x8a409296b6d78b70, 0xfd81a65ea4043ce1, 0x5dbf07d64599ff21, 0x797e5e0416255d02, 0x2da1f24e489f30a4}}, {{ 0xd1bbbce4d35711, 0xdba49d20e204088c, 0x6523385fd442067a, 0xfd08c31758f357c9, 0x45005c13ce4da178, 0xb447b9a98082cab8}, {0x60624f99df01e08, 0xa1b78d16406bbedb, 0x3ec8cd301a297e2b, 0xa3f13a8b99f0de0d, 0xa1f70d20526c3860, 0x141e5c2ab24a2cc6}}, {{ 0x62b2f9b766f8927, 0x84e89cf298f19686, 0x19b76882eb58d065, 0x83eb4a0974333090, 0xb1b258af14870a0f, 0x3a626e1d6c52a343}, {0xac00b71c66761d8, 0x65f002675c104942, 0x15312d3da70eb1ff, 0xe6ebc0389f7c74a9, 0x95db6bfd8d6c08b0, 0x1759cc53bd99e587}}, {{ 0xb2762e59e8b1726, 0x4f5f2cc230c408af, 0x946feeb304a05527, 0xac547a0fc23ecb0f, 0xd7b85f7fb50b7b2d, 0x330c869a5a6a172c}, {0xb9d02ed7832e13c, 0xc156e70185509cb0, 0x8832a5207a5b6073, 0x9b6752a6a9f063ab, 0x87f9cee6a6d13ab4, 0x3db1733518ed82ce}}, {{ 0x0fafbcfee42c127, 0x5365b3bb551050cb, 0x3e5af629dfaaf1f1, 0x2f7834c758e7963e, 0x240a7ebdaa847fb8, 0x1b4b4ea1b8e0a3f4}, {0xa3cc34905c12065, 0x2dfe831d67a0ca1f, 0x4ab3b00fbd1593b2, 0x648f561c2e2e9753, 0x9189bda376298faa, 0x31d6c32d20b7d611}}, {{ 0x145727f046231ba, 0x8beb1db86c744e42, 0x9b46e838effae6fb, 0x5fa2a9a03d5ab137, 0xd6fc56b57e6d388a, 0x3faf3d245796d7ca}, {0x4c62a0943bc6067, 0x5e16adb5aa0c5cd1, 0x7d982e8e6f8b8b4e, 0x1d91b64aacc66dd9, 0xa5590fe5de259cb8, 0x118433b551538046}}, {{ 0x0d3b237f4f996d, 0x44e7d73e1372f0ad, 0x9a3361efb8aec186, 0xbcea469bad92f5c9, 0x7cd44fe8b79b2c0a, 0x335c5cb38f901c0a}, {0x71d2cb72b6d6e7f, 0x83430219c88401d0, 0xa4c335245fd4a8ec, 0x07f85b9a9f2d4c93, 0xc152598aa9579dce, 0x3ffcca6bbbe3cb11}}, {{ 0xac4c88ea4bad10, 0x0987eb4f633da0f9, 0x2802264fff78d833, 0x52f0dff2c4af05e6, 0x577bdf7236ada3f6, 0x4cc5b2f06381c77c}, {0x9d0a6de9fd64e7b, 0x0b36d0975693a057, 0x4472326d12a1dee5, 0x0fd376c3bd2a7fef, 0xe3dc70f7de45eee2, 0x33d826828da72f84}}};

const digit_t p_cofactor_for_2f[NWORDS_P_COFACTOR_FOR_2F] = {0x3a8ffa868fe0799, 0x2daca7934d6edbe4, 0x824a11449f8c987f};
const uint16_t P_COFACTOR_FOR_2F_BITLENGTH = 188;
const uint16_t POWER_OF_2 = 191;
