# 前端部署指南 - GitHub Pages

本指南将帮助你将 Lotus Web Demo 前端部署到 GitHub Pages（免费）。

## 前置要求

1. GitHub 仓库已配置
2. 后端 API 已部署到 Render（参考 `DEPLOY_BACKEND.md`）

## 部署步骤

### 1. 更新 API 配置

在部署前端之前，确保后端已部署并获取了 Render URL。

1. 打开 `.github/workflows/web-demo.yml`
2. 找到这一行：
   ```yaml
   window.API_BASE_URL = window.API_BASE_URL || 'https://your-render-app.onrender.com';
   ```
3. 将 `your-render-app` 替换为你的实际 Render 应用名称
   - 例如：`https://lotus-api.onrender.com`

### 2. 启用 GitHub Pages

1. 进入你的 GitHub 仓库
2. 点击 **Settings** → **Pages**
3. 在 **Source** 下选择 **GitHub Actions**
4. 保存设置

### 3. 提交并推送代码

```bash
git add .
git commit -m "Configure GitHub Pages deployment"
git push origin main
```

### 4. 等待自动部署

1. 进入 GitHub 仓库的 **Actions** 标签
2. 查看 "Deploy Web Demo to GitHub Pages" workflow
3. 等待部署完成（通常 1-2 分钟）

### 5. 访问你的网站

部署完成后，你的网站将在以下 URL 可用：
```
https://<your-username>.github.io/Lotus/
```

或者如果你配置了自定义域名：
```
https://your-custom-domain.com
```

## 文件结构

部署后的文件结构：
- `index.html` - 主页面
- `app.js` - 应用逻辑
- `config.js` - API 配置（自动生成）

## 更新 API URL

如果需要更新后端 API URL：

1. 修改 `.github/workflows/web-demo.yml` 中的 API URL
2. 提交并推送更改
3. GitHub Actions 会自动重新部署

## 故障排除

### 部署失败

1. 检查 GitHub Actions 日志
2. 确保 `Lotus-Web-Demo` 目录存在且包含文件
3. 检查 workflow 文件语法是否正确

### API 连接失败

1. 检查 `config.js` 中的 API URL 是否正确
2. 确保后端服务正在运行（Render）
3. 检查浏览器控制台的错误信息
4. 确保后端 CORS 配置正确（已在代码中配置）

### 页面显示空白

1. 检查浏览器控制台是否有 JavaScript 错误
2. 检查网络请求是否成功
3. 确保所有文件都已正确部署

## 本地测试

在部署前，可以在本地测试：

1. 启动后端（本地）：
   ```bash
   python -m lotus.api.app
   ```

2. 在 `Lotus-Web-Demo` 目录创建 `config.js`：
   ```javascript
   window.API_BASE_URL = 'http://localhost:5000';
   ```

3. 在 `index.html` 中添加：
   ```html
   <script src="config.js"></script>
   ```

4. 使用本地服务器打开 `index.html`：
   ```bash
   cd Lotus-Web-Demo
   python -m http.server 8000
   ```

5. 访问 `http://localhost:8000`

## 下一步

- 配置自定义域名（可选）
- 设置 HTTPS（GitHub Pages 自动提供）
- 监控使用情况

