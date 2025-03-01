// Generate a markdown file with a table summary of the aggregated benchmarks.
// If a benchmark-base file is available, shows the comparison against base (ie master in a PR).
import { createConsoleLogger } from '@aztec/foundation/log';
import { BENCHMARK_HISTORY_BLOCK_SIZE, Metrics } from '@aztec/types/stats';

import * as fs from 'fs';
import pick from 'lodash.pick';

import { BaseBenchFile, BenchFile } from './paths.js';

// Input file paths
const inputFile = BenchFile;
const baseFile = BaseBenchFile;

const COMMENT_MARK = '<!-- AUTOGENERATED BENCHMARK COMMENT -->';
const S3_URL = 'https://aztec-ci-artifacts.s3.us-east-2.amazonaws.com';

const log = createConsoleLogger();

/** Returns a cell content formatted as string */
function getCell(
  data: Record<string, Record<string, number>>,
  base: Record<string, Record<string, number>> | undefined,
  row: string,
  col: string,
) {
  const value = data[row][col];
  const baseValue = base ? (base[row] ?? {})[col] : undefined;
  const percentDiff = baseValue ? Math.round(((value - baseValue) / baseValue) * 100) : undefined;
  const formattedValue = formatValue(value);
  const highlight = percentDiff && Math.abs(percentDiff) > 10 ? '**' : '';
  const warning = percentDiff && Math.abs(percentDiff) > 10 ? ':warning:' : '';
  const percentSign = percentDiff && percentDiff > 0 ? '+' : '';
  return percentDiff && Math.abs(percentDiff) >= 1
    ? `${warning} ${formattedValue} (${highlight}<span title="${formatValue(
        baseValue!,
      )}">${percentSign}${percentDiff}%</span>${highlight})`
    : formattedValue;
}

/** Returns the description of a metric name, if found. */
function tryGetDescription(name: string) {
  return Metrics.find(m => m.name === name)?.description;
}

/** Wraps the metric name in a span with a title with the description, if found. */
function withDescriptionTitle(name: string) {
  const description = tryGetDescription(name);
  if (!description) return name;
  return `<span title="${description}">${name}</span>`;
}

/** Formats a numeric value for display. */
function formatValue(value: number) {
  if (value < 100) return value.toPrecision(3);
  return value.toLocaleString();
}

/** Transposes an object topmost and nested keys. */
function transpose(obj: any) {
  const transposed: any = {};
  for (const outerKey in obj) {
    const innerObj = obj[outerKey];
    for (const innerKey in innerObj) {
      if (!transposed[innerKey]) transposed[innerKey] = {};
      transposed[innerKey][outerKey] = innerObj[innerKey];
    }
  }
  return transposed;
}

/** Returns the base benchmark for comparison, if exists */
function getBaseBenchmark(): Record<string, Record<string, number>> | undefined {
  try {
    return JSON.parse(fs.readFileSync(baseFile, 'utf-8'));
  } catch {
    return undefined;
  }
}

/** Creates a table in md out of the data (rows and cols). */
function getTableContent(
  data: Record<string, Record<string, number>>,
  baseBenchmark: Record<string, Record<string, number>> | undefined,
  groupUnit = '',
  col1Title = 'Metric',
) {
  const rowKeys = Object.keys(data);
  const groups = [...new Set(rowKeys.flatMap(key => Object.keys(data[key])))];
  const makeHeader = (colTitle: string) => `${withDescriptionTitle(colTitle)} ${groupUnit}`;
  const header = `| ${col1Title} | ${groups.map(makeHeader).join(' | ')} |`;
  const separator = `| - | ${groups.map(() => '-').join(' | ')} |`;
  const makeCell = (row: string, col: string) => getCell(data, baseBenchmark, row, col);
  const rows = rowKeys.map(key => `${withDescriptionTitle(key)} | ${groups.map(g => makeCell(key, g)).join(' | ')} |`);

  return `
${header}
${separator}
${rows.join('\n')}
  `;
}

/** Creates a md with the benchmark contents. */
export function getMarkdown() {
  const benchmark = JSON.parse(fs.readFileSync(inputFile, 'utf-8'));
  const baseBenchmark = getBaseBenchmark();
  const metricsByBlockSize = Metrics.filter(m => m.groupBy === 'block-size').map(m => m.name);
  const metricsByChainLength = Metrics.filter(m => m.groupBy === 'chain-length').map(m => m.name);
  const metricsByCircuitName = Metrics.filter(m => m.groupBy === 'circuit-name').map(m => m.name);

  const baseHash = process.env.BASE_COMMIT_HASH;
  const baseUrl = baseHash && `[\`${baseHash.slice(0, 8)}\`](${S3_URL}/benchmarks-v1/master/${baseHash}.json)`;
  const baseCommitText = baseUrl
    ? `\nValues are compared against data from master at commit ${baseUrl} and shown if the difference exceeds 1%.`
    : '';

  const prNumber = process.env.CIRCLE_PULL_REQUEST && parseInt(process.env.CIRCLE_PULL_REQUEST.split('/')[6]);
  const prSourceDataUrl = prNumber && `${S3_URL}/benchmarks-v1/pulls/${prNumber}.json`;
  const prSourceDataText = prSourceDataUrl
    ? `\nThis benchmark source data is available in JSON format on S3 [here](${prSourceDataUrl}).`
    : '';

  return `
## Benchmark results

All benchmarks are run on txs on the \`Benchmarking\` contract on the repository. Each tx consists of a batch call  to \`create_note\` and \`increment_balance\`, which guarantees that each tx has a private call, a nested private call, a public call, and a nested public call, as well as an emitted private note, an unencrypted log, and public storage read and write. 
${prSourceDataText}
${baseCommitText}

### L2 block published to L1

Each column represents the number of txs on an L2 block published to L1.
${getTableContent(pick(benchmark, metricsByBlockSize), baseBenchmark, 'txs')}

### L2 chain processing

Each column represents the number of blocks on the L2 chain where each block has ${BENCHMARK_HISTORY_BLOCK_SIZE} txs.
${getTableContent(pick(benchmark, metricsByChainLength), baseBenchmark, 'blocks')}

### Circuits stats

Stats on running time and I/O sizes collected for every circuit run across all benchmarks.
${getTableContent(transpose(pick(benchmark, metricsByCircuitName)), transpose(baseBenchmark), '', 'Circuit')}

${COMMENT_MARK}
`;
}

/** Entrypoint */
export function main() {
  log(getMarkdown());
}
