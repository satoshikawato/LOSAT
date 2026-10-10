// The app's version (package.json) and the git commit it was built from (short SHA, or
// "unknown"), set by vite.config.ts and written into session files. A build that does not
// define them (the E2E harness's) leaves them undefined; read them with `typeof`.
declare const __LOSAT_APP_VERSION__: string;
declare const __LOSAT_APP_BUILD__: string;
