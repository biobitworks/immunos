# IMMUNOS Dashboard Feature Roadmap

**Last Reviewed**: 2025-12-25  
**Current Version**: v20251225-kb8  
**Status**: Phase 2 (Training Expansion) - Active

---

## Executive Summary

The IMMUNOS Dashboard is a comprehensive web interface for monitoring and managing the IMMUNOS multi-agent system. This document identifies prudent improvements, security enhancements, and feature additions based on current implementation review.

---

## 🔴 Critical Security Issues (High Priority)

### 1. Hardcoded Secret Key
**Current**: `SECRET_KEY = 'immunos-dashboard-secret-key-change-in-production'`  
**Risk**: High - Session hijacking, CSRF attacks  
**Fix**:
- Load from environment variable: `os.getenv('IMMUNOS_SECRET_KEY')`
- Generate random key on first run if not set
- Store in `.immunos/config/dashboard_secret.key` (gitignored)

**Priority**: 🔴 **IMMEDIATE**

### 2. CORS Wide Open
**Current**: `CORS(app)` and `cors_allowed_origins="*"`  
**Risk**: Medium - Cross-origin attacks  
**Fix**:
- Environment-based CORS config
- Default to localhost only in production
- Configurable allowed origins via config file

**Priority**: 🔴 **IMMEDIATE**

### 3. CDN Dependencies (Air-Gapped)
**Current**: Tailwind, Chart.js, Socket.IO from CDN  
**Risk**: Medium - Won't work offline/air-gapped  
**Fix**:
- Vendor all JS/CSS locally
- Create `static/vendor/` directory
- Download and include in bundle

**Priority**: 🟡 **HIGH** (for air-gapped deployment)

---

## 🟡 Important Features (Medium Priority)

### 4. MCP Server Integration
**Current**: Dashboard doesn't integrate with IMMUNOS-MCP server  
**Opportunity**: 
- Add MCP server status indicator
- Show active MCP tools (immune_scan, immune_anomaly, etc.)
- Display orchestrator mode (local/remote)
- Show MCP connection health

**Priority**: 🟡 **HIGH** (aligns with new MCP architecture)

### 5. Local Thymus/Orchestrator Status
**Current**: No visibility into local orchestrator state  
**Feature**:
- Display current orchestrator mode (air-gapped/online)
- Show fallback status (local vs remote)
- Rate limit detection indicator
- Network connectivity status

**Priority**: 🟡 **HIGH** (core functionality)

### 6. Training Result Parsing
**Current**: Training outputs not parsed into history table  
**Status**: Partially complete (Phase 2)  
**Feature**:
- Parse training outputs → populate detectors/accuracy
- Display training metrics in Thymus tab
- Show accuracy trends over time
- Export training reports

**Priority**: 🟡 **MEDIUM** (completes Phase 2)

### 7. Issue Auto-Tagging
**Current**: Manual issue creation  
**Feature**:
- Auto-create issues from detection events
- Tag issues based on anomaly patterns
- Link issues to specific code/files
- Priority assignment based on confidence

**Priority**: 🟡 **MEDIUM** (Phase 3)

---

## 🟢 Enhancement Features (Lower Priority)

### 8. Health Check Endpoints
**Current**: No comprehensive health checks  
**Feature**:
- `/api/health` endpoint with component status
- Ollama connectivity check
- Database writeability check
- Agent artifact availability
- Model availability status

**Priority**: 🟢 **MEDIUM**

### 9. Authentication/Authorization
**Current**: No auth system  
**Feature**:
- Optional authentication for multi-user setups
- Shared token gating for `/api/*` endpoints
- Role-based access control
- Session management

**Priority**: 🟢 **LOW** (only needed for shared machines)

### 10. Daily/Weekly Summary Export
**Current**: No automated summaries  
**Feature**:
- Daily summary export to `daily/YYYY-MM-DD.md`
- Weekly aggregation reports
- Anomaly detection summaries
- Training progress reports

**Priority**: 🟢 **MEDIUM** (Phase 3)

### 11. Standardized Logging Schema
**Current**: Inconsistent logging across components  
**Feature**:
- Unified logging schema
- Structured logging (JSON)
- Log aggregation endpoint
- Search/filter capabilities

**Priority**: 🟢 **MEDIUM** (Phase 3)

### 12. Cost/Latency Based Routing
**Current**: Manual orchestrator selection  
**Feature**:
- Automatic cost calculation
- Latency-based routing
- Backend profile management
- Routing strategy configuration

**Priority**: 🟢 **LOW** (Phase 4)

### 13. Offline Install Script
**Current**: Manual installation  
**Feature**:
- Automated offline installation script
- Ubuntu USB image support
- Dependency bundling
- Configuration wizard

**Priority**: 🟢 **LOW** (Phase 5)

### 14. Lite vs Full Runtime Profiles
**Current**: Single configuration  
**Feature**:
- Lite profile (minimal features, low resource)
- Full profile (all features, full resources)
- Profile switching
- Resource monitoring

**Priority**: 🟢 **LOW** (Phase 5)

---

## 📊 Dashboard-Specific Enhancements

### 15. Real-Time MCP Tool Status
**Feature**: Show which MCP tools are active/available
- Tool availability indicator
- Last used timestamp
- Success/failure rates
- Tool-specific metrics

**Priority**: 🟡 **HIGH**

### 16. Orchestrator Mode Toggle UI
**Feature**: Visual toggle for air-gapped/online mode
- Mode switch button
- Current mode indicator
- Fallback status display
- Network connectivity indicator

**Priority**: 🟡 **HIGH**

### 17. Agent Health Dashboard
**Feature**: Comprehensive agent status view
- B Cell pattern matching stats
- NK Cell anomaly detection stats
- Dendritic feature extraction stats
- Memory agent cache stats
- Per-agent performance metrics

**Priority**: 🟡 **MEDIUM**

### 18. Detection Event Timeline
**Feature**: Visual timeline of detection events
- Chronological event display
- Filter by agent type
- Filter by severity
- Export event logs

**Priority**: 🟢 **MEDIUM**

### 19. Configuration Management UI
**Feature**: Web-based config editor
- Edit `~/.immunos-mcp/config.yaml` from UI
- Validate configuration
- Test configuration changes
- Import/export configs

**Priority**: 🟢 **LOW**

### 20. Model Management Integration
**Feature**: Better Ollama model management
- Model download progress
- Model size/performance info
- Model usage statistics
- Model recommendation engine

**Priority**: 🟢 **LOW**

---

## 🔧 Technical Improvements

### 21. API Documentation
**Current**: No API docs  
**Feature**:
- OpenAPI/Swagger spec
- Interactive API explorer
- Endpoint documentation
- Request/response examples

**Priority**: 🟢 **MEDIUM**

### 22. Error Handling & Logging
**Current**: Basic error handling  
**Feature**:
- Comprehensive error logging
- Error notification system
- Error recovery mechanisms
- User-friendly error messages

**Priority**: 🟡 **MEDIUM**

### 23. Database Migration System
**Current**: Manual schema updates  
**Feature**:
- Versioned migrations
- Automatic migration on startup
- Migration rollback capability
- Migration status tracking

**Priority**: 🟢 **LOW**

### 24. Performance Monitoring
**Current**: No performance metrics  
**Feature**:
- Response time tracking
- Database query performance
- WebSocket connection metrics
- Resource usage monitoring

**Priority**: 🟢 **LOW**

### 25. Testing Infrastructure
**Current**: No automated tests  
**Feature**:
- Unit tests for API endpoints
- Integration tests for dashboard
- E2E tests for critical flows
- Smoke tests for deployment

**Priority**: 🟡 **MEDIUM** (Phase 6)

---

## 📱 User Experience Improvements

### 26. Responsive Design
**Current**: Desktop-focused  
**Feature**:
- Mobile-responsive layout
- Tablet optimization
- Touch-friendly controls
- Adaptive UI components

**Priority**: 🟢 **LOW**

### 27. Dark/Light Theme Toggle
**Current**: Dark theme only  
**Feature**:
- Light theme option
- Theme persistence
- System theme detection
- Custom theme support

**Priority**: 🟢 **LOW**

### 28. Keyboard Shortcuts
**Current**: Mouse-only navigation  
**Feature**:
- Keyboard navigation
- Shortcut key bindings
- Command palette (Cmd+K)
- Quick actions

**Priority**: 🟢 **LOW**

### 29. Export/Import Functionality
**Current**: Limited export options  
**Feature**:
- Export detection results (CSV/JSON)
- Export configuration
- Import detection data
- Bulk operations

**Priority**: 🟢 **MEDIUM**

### 30. Search & Filtering
**Current**: Basic filtering  
**Feature**:
- Global search
- Advanced filters
- Saved filter presets
- Search history

**Priority**: 🟢 **MEDIUM**

---

## 🚀 Integration Features

### 31. CI/CD Integration
**Feature**: Webhook support for CI/CD
- GitHub Actions integration
- GitLab CI integration
- Automated scanning on commits
- PR comment integration

**Priority**: 🟢 **LOW**

### 32. IDE Plugin Support
**Feature**: IDE integration
- VS Code extension
- Cursor integration (already via MCP)
- JetBrains plugin
- Command-line interface

**Priority**: 🟢 **LOW**

### 33. Notification System
**Feature**: Alert notifications
- Email notifications
- Desktop notifications
- Webhook notifications
- Slack/Discord integration

**Priority**: 🟢 **LOW**

---

## 📋 Implementation Priority Matrix

### Immediate (This Week)
1. ✅ Fix hardcoded SECRET_KEY
2. ✅ Fix CORS configuration
3. ✅ Add MCP server status indicator
4. ✅ Add local orchestrator status display

### Short-Term (This Month)
5. ✅ Vendor CDN dependencies
6. ✅ Add health check endpoints
7. ✅ Implement training result parsing
8. ✅ Add orchestrator mode toggle UI

### Medium-Term (Next Quarter)
9. ✅ Issue auto-tagging
10. ✅ Daily/weekly summary export
11. ✅ Standardized logging schema
12. ✅ Agent health dashboard

### Long-Term (Future)
13. ✅ Authentication system
14. ✅ Cost/latency routing
15. ✅ Offline install script
16. ✅ Lite/Full profiles

---

## 🎯 Feature Categories

### Security & Compliance
- Hardcoded secrets fix
- CORS configuration
- Authentication (optional)
- Audit logging

### Core Functionality
- MCP server integration
- Local orchestrator status
- Training result parsing
- Issue auto-tagging

### User Experience
- Health checks
- Configuration UI
- Export/import
- Search & filtering

### Technical Debt
- CDN dependencies
- API documentation
- Testing infrastructure
- Database migrations

### Integration
- CI/CD webhooks
- IDE plugins
- Notification system
- External tool integration

---

## 📝 Notes

- **Current Phase**: Phase 2 (Training Expansion)
- **Next Milestone**: Complete training result parsing
- **Blockers**: None identified
- **Dependencies**: IMMUNOS-MCP server (completed)

---

## 🔄 Review Schedule

- **Weekly**: Review critical security issues
- **Monthly**: Review feature priorities
- **Quarterly**: Major roadmap update

---

**Last Updated**: 2025-12-25  
**Next Review**: 2026-01-01

