"""Deterministic network-only httpx stub for tests, with no optional respx dependency.

Intercepts only outbound AsyncClient instances created without an explicit transport.
ASGI in-process test clients with their own transport are unaffected.
Unexpected external HTTP calls fail closed instead of leaving the test machine.
"""
from __future__ import annotations

from collections.abc import Callable
from urllib.parse import urlsplit

import httpx


class Route:
    def __init__(self):
        self.response: httpx.Response | None = None
        self.callback: Callable | None = None

    def mock(self, *, return_value: httpx.Response | None = None, side_effect: Callable | None = None):
        if (return_value is None) == (side_effect is None):
            raise ValueError("Provide exactly one HTTP response or callback")
        self.response, self.callback = return_value, side_effect
        return self

    def execute(self, request: httpx.Request) -> httpx.Response:
        if self.callback is not None:
            return self.callback(request)
        if self.response is not None:
            # Create an independent response object so repeated calls do not
            # reuse a consumed stream or an old request context.
            return httpx.Response(
                self.response.status_code,
                headers=dict(self.response.headers),
                content=self.response.content,
            )
        raise AssertionError(f"Unconfigured mock request: {request.method} {request.url}")


class Router:
    def __init__(self, origin: str):
        self.origin = origin.rstrip("/")
        parsed = urlsplit(self.origin)
        self.scheme = parsed.scheme
        self.netloc = parsed.netloc
        self.routes: dict[tuple[str, str], Route] = {}
        self.called: list[tuple[str, str]] = []

    def register(self, method: str, path: str) -> Route:
        key = (method, path)
        return self.routes.setdefault(key, Route())

    def get(self, path: str) -> Route:
        return self.register("GET", path)

    def post(self, path: str) -> Route:
        return self.register("POST", path)

    def respond(self, request: httpx.Request) -> httpx.Response:
        url = request.url
        if url.scheme != self.scheme or url.netloc.decode() != self.netloc:
            raise AssertionError(f"Unexpected network host: {request.method} {request.url}")
        key = (request.method, url.path)
        self.called.append(key)
        route = self.routes.get(key)
        if route is None:
            raise AssertionError(f"Unexpected external HTTP request: {request.method} {request.url}")
        return route.execute(request)


def install(monkeypatch, origin: str) -> Router:
    """Use httpx.MockTransport for new outbound AsyncClients only."""
    router = Router(origin)
    previous = httpx.AsyncClient.__init__

    def patched(self, *args, **kwargs):
        if kwargs.get("transport") is None:
            kwargs["transport"] = httpx.MockTransport(router.respond)
        return previous(self, *args, **kwargs)

    monkeypatch.setattr(httpx.AsyncClient, "__init__", patched)
    return router
