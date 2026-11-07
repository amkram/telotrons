// Build SUPPLEMENTARY_TABLES.docx and SUPPLEMENTARY_FIGURES.docx under $TELO_PACKAGE/supplementary
// Usage: node build_supplement_docx.js [tables|figures], both without an argument.
// Needs `npm install docx`.
const fs = require('fs');
const path = require('path');
const {
  Document,
  Packer,
  Paragraph,
  TextRun,
  Table,
  TableRow,
  TableCell,
  WidthType,
  ShadingType,
  PageOrientation,
  AlignmentType,
  PageBreak,
  BorderStyle,
  ImageRun
} = require('docx');

const PACKAGE = process.env.TELO_PACKAGE;
const TABLE_DIR = path.join(PACKAGE, 'supplementary', 'tables');
const FIGURE_DIR = path.join(PACKAGE, 'figures');
const which = process.argv[2] || 'both';

// rows of a TSV under supplementary/tables, split into fields
const tsv = f =>
  fs.readFileSync(path.join(TABLE_DIR, f), 'utf8').trim().split('\n').map(l => l.split('\t'));

// table width in twips: landscape letter minus half-inch margins
const TABLE_WIDTH = 14400;
const border = {style: BorderStyle.SINGLE, size: 4, color: '808080'};
const borders = {top: border, bottom: border, left: border, right: border};

// one table cell holding a single Arial run (size in half-points)
function cell(text, w, {bold = false, size = 14, shade = null, align = AlignmentType.LEFT} = {}) {
  return new TableCell({
    width: {size: w, type: WidthType.DXA},
    borders,
    margins: {top: 30, bottom: 30, left: 50, right: 50},
    shading: shade ? {type: ShadingType.CLEAR, fill: shade, color: 'auto'} : undefined,
    children: [
      new Paragraph({
        alignment: align,
        spacing: {before: 0, after: 0},
        children: [new TextRun({text, bold, size, font: 'Arial'})]
      })
    ]
  });
}

// a table with a repeated, shaded header row, its column widths adding up to TABLE_WIDTH
function table(rows, widths, {size = 14, numericFrom = null, totalRow = false} = {}) {
  const sum = widths.reduce((a, b) => a + b, 0);
  if (sum !== TABLE_WIDTH) throw new Error('widths sum ' + sum);
  if (rows[0].length !== widths.length) {
    throw new Error(
      'column count ' + rows[0].length + ' vs widths ' + widths.length +
        ' for header ' + rows[0].join('|')
    );
  }

  const [headerCells, ...body] = rows;
  const header = new TableRow({
    tableHeader: true,
    cantSplit: true,
    children: headerCells.map((t, i) => cell(t, widths[i], {bold: true, size, shade: 'E7E6E6'}))
  });

  // columns from numericFrom on are right-aligned, and totalRow bolds the last row
  const bodyRows = body.map((r, ri) => new TableRow({
    cantSplit: true,
    children: r.map((t, i) => cell(t, widths[i], {
      size,
      bold: totalRow && ri === body.length - 1,
      align: (numericFrom !== null && i >= numericFrom) ? AlignmentType.RIGHT : AlignmentType.LEFT
    }))
  }));

  return new Table({
    width: {size: TABLE_WIDTH, type: WidthType.DXA},
    columnWidths: widths,
    rows: [header, ...bodyRows]
  });
}

// paragraph helpers: body text, table title, part heading, page break
const para = (t, {size = 18, bold = false, before = 120, after = 120} = {}) => new Paragraph({
  spacing: {before, after},
  children: [new TextRun({text: t, size, bold, font: 'Arial'})]
});
const title = t => para(t, {size: 20, bold: true, before: 0, after: 120});
const heading = t => para(t, {size: 18, bold: true, before: 240, after: 120});
const pageBreak = () => new Paragraph({children: [new PageBreak()]});

// a landscape letter section with half-inch margins
const landscape = children => ({
  properties: {
    page: {
      size: {width: 12240, height: 15840, orientation: PageOrientation.LANDSCAPE},
      margin: {top: 720, bottom: 720, left: 720, right: 720}
    }
  },
  children
});

// table captions
const CAPTIONS = {
  S1:
    'Telotron-bearing genomes and their controls. Part A identifies each genome: source, whether ' +
    'introns were taken from annotated exons or CDS blocks, the role in the comparison, the ' +
    'telomere motif of its lineage, the number of annotated introns screened, the telotron ' +
    'counts at the two coverage thresholds and, for MAGs, BUSCO completeness from the SMAG ' +
    'release (Delmont et al. 2022). Every genome was screened with the published telomere motif ' +
    'of its lineage in all cyclic rotations and both orientations (TTAGGG for haptophytes, the ' +
    'chrysophyte and the other stramenopile; TTTAGGG for dinoflagellates, Eimeria and the ' +
    'cryptophyte; TTAGGC for the nematode); controls have no intron above 50% repeat coverage ' +
    'under that motif. Part B gives the per-genome results: primary loci whose closest ' +
    'qualifying ortholog lacks the intron over the loci with a qualifying comparison (intron ' +
    'absent or an annotation-supported ordinary intron at the same position and phase), putative ' +
    'target-site duplications at primary loci, and interstitial telomeric sequence (ITS) density ' +
    'outside genes with the eligible sequence and the paired genome used in Figure 3.',
  S2:
    'The complete catalogue of 51,431 telotrons is provided as a separate spreadsheet ' +
    '(SUPPLEMENTARY_TABLES.xlsx, sheet S2_catalogue). Each row gives coordinates, strand, intron ' +
    'length, splice boundaries, repeat coverage (primary loci are those above 0.85), ' +
    'architecture and the first rule failed by unclassified loci, linker coordinates, putative ' +
    'TSD length and sequence, RNA junction reads where RNA-seq was available, the closest ' +
    'qualifying comparator state with the number of eligible comparators, and the exact position ' +
    'and phase audit of ordinary-intron comparisons.',
  S3:
    'RNA-seq reads were aligned to the genome assembly with minimap2 (splice:sr preset, ' +
    'secondary alignments suppressed). A locus counts as supported when at least one primary ' +
    'alignment with MAPQ ≥20 spans the exact annotated junction with ≥10 contiguous aligned ' +
    'bases on each side; counts are alignment records, not independent molecules. The last ' +
    'column gives the percentage supported among a fixed random sample of 2,000 ordinary ' +
    'annotated introns per genome (≥30 bp, below 30% telomeric coverage, not in the catalogue). ' +
    'Library types are as registered in ENA; SRR485952 is an EST library.'
};

// tables document: S1 parts A and B, the S2 note, S3
if (which !== 'figures') {
  const s1PartA = tsv('S1_genomes_A.tsv');
  const s1PartB = tsv('S1_genomes_B.tsv');
  const s3 = tsv('S3_rna.tsv');

  // column widths in twips
  const tables = new Document({
    styles: {default: {document: {run: {font: 'Arial', size: 18}}}},
    sections: [landscape([
      title('Supplementary Table S1. Genomes: identity, screen status and per-genome results.'),
      para(CAPTIONS.S1),
      heading('A. Genomes and screens'),
      table(
        s1PartA,
        [1450, 3000, 2300, 1250, 950, 800, 950, 1150, 800, 800, 950],
        {size: 12, numericFrom: 7}
      ),
      pageBreak(),
      heading('B. Per-genome results'),
      table(s1PartB, [1500, 3200, 2500, 1900, 1200, 1300, 1200, 1600], {size: 12, numericFrom: 3}),
      pageBreak(),
      title('Supplementary Table S2. Telotron catalogue with per-locus evidence.'),
      para(CAPTIONS.S2),
      title(
        'Supplementary Table S3. RNA-seq datasets, junction support, ' +
          'and the ordinary-intron baseline.'
      ),
      para(CAPTIONS.S3),
      table(
        s3,
        [1900, 1500, 2900, 1700, 1100, 1100, 1200, 1200, 1800],
        {size: 12, numericFrom: 4, totalRow: true}
      )
    ])]
  });

  Packer.toBuffer(tables).then(buffer => {
    fs.writeFileSync(path.join(PACKAGE, 'supplementary', 'SUPPLEMENTARY_TABLES.docx'), buffer);
    console.log('tables docx', buffer.length);
  });
}

// figures document: one portrait page per supplementary figure
if (which !== 'tables') {
  // caption text from figures/<stem>.caption.md, without bold markers
  const figureCaption = stem =>
    fs.readFileSync(path.join(FIGURE_DIR, stem + '.caption.md'), 'utf8')
      .trim()
      .replace(/\*\*/g, '');

  // a centred PNG, widthInches wide at 96 px per inch, height from its aspect ratio
  // (pixel width and height are bytes 16-23 of the PNG header)
  const image = (stem, widthInches) => {
    const data = fs.readFileSync(path.join(FIGURE_DIR, stem + '.png'));
    const heightInches = widthInches * data.readUInt32BE(20) / data.readUInt32BE(16);
    return new Paragraph({
      alignment: AlignmentType.CENTER,
      children: [new ImageRun({
        type: 'png',
        data,
        transformation: {width: Math.round(widthInches * 96), height: Math.round(heightInches * 96)}
      })]
    });
  };

  // caption paragraph with its first sentence in bold
  const captionParagraph = t => {
    const i = t.indexOf('. ');
    const head = t.slice(0, i + 1);
    const rest = t.slice(i + 1).trim();
    return new Paragraph({
      spacing: {before: 160, after: 120},
      children: [
        new TextRun({text: head + ' ', bold: true, size: 18, font: 'Arial'}),
        new TextRun({text: rest, size: 18, font: 'Arial'})
      ]
    });
  };

  // a portrait letter section with half-inch margins
  const portrait = children => ({
    properties: {
      page: {
        size: {width: 12240, height: 15840},
        margin: {top: 720, bottom: 720, left: 720, right: 720}
      }
    },
    children
  });

  // figure stems and printed widths in inches
  const figures = new Document({
    styles: {default: {document: {run: {font: 'Arial', size: 18}}}},
    sections: [
      ['FigureS1_unclassified_reasons', 7.0],
      ['FigureS2_genome_tree', 7.0],
      ['FigureS3_single_motif_introns', 7.0]
    ].map(([stem, width]) => portrait([image(stem, width), captionParagraph(figureCaption(stem))]))
  });

  Packer.toBuffer(figures).then(buffer => {
    fs.writeFileSync(path.join(PACKAGE, 'supplementary', 'SUPPLEMENTARY_FIGURES.docx'), buffer);
    console.log('figures docx', buffer.length);
  });
}
