# Hash scope notes (2026-10-08 audit)

This sidecar note labels what the hash and checksum fields in the files below actually cover. The data files themselves are unchanged. It is a documentation clarification, not a license verdict and not a source re-admission.

Audited commit: `fa557b0ee057511bd459527cc2f4c7fb6c47d648`

## `results/dna_fountain_source_audit.json`
- Lines (verified against the audited commit): 8, 13, 19, 34, 50, 66, 82, 98, 114, 130, 146, 162, 178, 194, 210, 226, 242, 258, 274 (19 lines)
- Scope: The tsv_sha256 values hash ENA report/metadata tables, and any md5 is publisher-advertised. No FASTQ payload was acquired, so these are not payload integrity checks.

## `data/payloads/manifest.json`
- Lines (verified against the audited commit): 3, 5, 8, 10, 13, 15, 18, 20, 23, 25, 28, 30, 33, 35, 38, 40, 43, 45, 48, 50, 53, 55, 58, 60, 63, 65, 68, 70, 73, 75, 78, 80, 83, 85, 88, 90, 93, 95, 98, 100, 103, 105, 108, 110, 113, 115, 118, 120, 123, 125, 128, 130, 133, 135, 138, 140, 143, 145, 148, 150, 153, 155, 158, 160, 163, 165, 168, 170, 173, 175, 178, 180, 183, 185, 188, 190, 193, 195, 198, 200, 203, 205, 208, 210, 213, 215, 218, 220, 223, 225, 228, 230, 233, 235, 238, 240, 243, 245, 248, 250, 253, 255, 258, 260, 263, 265, 268, 270, 273, 275, 278, 280, 283, 285, 288, 290, 293, 295, 298, 300, 303, 305, 308, 310, 313, 315, 318, 320, 323, 325, 328, 330, 333, 335, 338, 340, 343, 345, 348, 350, 353, 355, 358, 360, 363, 365, 368, 370, 373, 375, 378, 380, 383, 385, 388, 390, 393, 395, 398, 400, 403, 405, 408, 410, 413, 415, 418, 420, 423, 425, 428, 430, 433, 435, 438, 440, 443, 445, 448, 450, 453, 455, 458, 460, 463, 465, 468, 470, 473, 475, 478, 480, 483, 485, 488, 490, 493, 495, 498, 500, 503, 505, 508, 510, 513, 515, 518, 520, 523, 525, 528, 530, 533, 535, 538, 540, 543, 545, 548, 550, 553, 555, 558, 560, 563, 565, 568, 570, 573, 575, 578, 580, 583, 585, 588, 590, 593, 595, 598, 600, 603, 605, 608, 610, 613, 615, 618, 620, 623, 625, 628, 630, 633, 635, 638, 640, 643, 645, 648, 650 (260 lines)
- Scope: This manifest lists accessions and lengths only. It carries no source-payload hash.

