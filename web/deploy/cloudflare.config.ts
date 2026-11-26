// Cloudflare Workers Static Assets configuration for LOSAT Web (`cf` CLI).
//
// Built static files in `web/app/dist/` (including `public/_headers`, `/engine/*.wasm`, and
// license texts) are deployed as an assets-only Worker with:
//   - Preview:    `cf previews deploy ../app/dist --preview-name <name>`
//   - Production: `cf deploy ../app/dist`
//
// Custom domains (such as `losat.app`) are attached only on non-preview deploys when
// `LOSAT_WEB_DOMAINS` is set (comma-separated hostnames), so preview deploys never depend on
// DNS zone onboarding. Preview URLs stay enabled (`previewUrls: true`) and Cloudflare sends
// `X-Robots-Tag: noindex` automatically on `workers.dev` preview URLs.

interface ConfigContext {
  readonly isPreview?: boolean;
  readonly mode?: string;
}

function parseDomains(raw: string | undefined): string[] {
  if (!raw) return [];
  return raw
    .split(',')
    .map((entry) => entry.trim())
    .filter((entry) => entry.length > 0);
}

export default ({ isPreview = false }: ConfigContext = {}) => {
  const workerName = (process.env.LOSAT_WEB_WORKER_NAME ?? 'losat-web').trim() || 'losat-web';
  const domains = isPreview ? [] : parseDomains(process.env.LOSAT_WEB_DOMAINS);
  const workersDev = isPreview ? true : domains.length === 0 || process.env.LOSAT_WEB_WORKERS_DEV === 'true';

  return {
    worker: {
      name: workerName,
      compatibilityDate: '2026-11-26',
      previewUrls: true,
      workersDev,
      ...(domains.length > 0 ? { domains } : {}),
      assets: {
        htmlHandling: 'auto-trailing-slash',
        notFoundHandling: 'none',
      },
    },
  };
};
