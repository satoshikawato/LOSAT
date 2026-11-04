- Web tests: every navigation of the E2E tests (`page.goto`, `page.reload`) goes through one helper that
  waits 15 s and then navigates once more, recording a `navigation made again` annotation. WebKit now and
  then never received a page's document response (21 `page.goto` timeouts in CI in 7 days, 4 in 4,240
  local navigations); the second navigation got it within a second each time. The application is unchanged.
